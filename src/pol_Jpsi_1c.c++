#include "kinematics/RamboEP.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>

#include <tclap/CmdLine.h>


dbl_type Mcc = 3.1, mc = Mcc / 2;
dbl_type O3S11 = 0.270;  // from Gee 
dbl_type O3S18 = 6.6e-3; // from hep-ph/9511315, Table I
dbl_type O3P08 = 2.2e-2*mc*mc; // from hep-ph/9511315, Table II
dbl_type O1S08 = 2.2e-2;
dbl_type Opsi=0.270;
dbl_type nanob = 0.389e6, picob = 1e3 * nanob; // conversion to barn
dbl_type PI = acos(-1), alpha = 1. / 137, alphas = 0.3;
dbl_type x;



#include "LHAPDF/LHAPDF.h"
LHAPDF::PDF *lhapdf_pdf;

void write_histogram_to_file(TH1F &histogram, string file_name);
bool kinematics(RamboEP *ramEP);
dbl_type getMatr2_3S1_cs();
dbl_type getMatr2_3S1_co();
dbl_type getMatr2_1S0_co();
dbl_type getMatr2_3P0_co();

dbl_type P[4], k[4], kp[4], k1[4], k2[4], k3[4], pPsi[4];

int nEv;
dbl_type ecm;
dbl_type Q2, minQ2, maxQ2;
dbl_type W2, minW2, maxW2;
dbl_type z, zMin, zMax;
dbl_type Y;
bool gauge;
string pdf_set_name;
int nBins;
dbl_type xi;
bool ldme;
string out_prefix;
bool saveHist;
int seed;
bool check;




void init_from_command_line(int argc, char **argv) {
    cout << "Started with command" << endl;
    for (int i = 0; i < argc; ++i) cout << argv[i] << " ";
    cout << endl;

    TCLAP::CmdLine cmd("J/psi meson photoproduction MC simulator", ' ', "1.0");
    TCLAP::ValueArg<double> arg_nEv("n", "nEv", "number of events to generate", false, 1e7, "double", cmd);
    TCLAP::ValueArg<double> arg_ecm("e", "ecm", "reaction energy in cms frame [GeV]", false, 300., "double", cmd);
    TCLAP::ValueArg<double> arg_minQ2("", "minQ2", "minimal Q2 value [GeV^2]", false, 4., "double", cmd);
    TCLAP::ValueArg<double> arg_maxQ2("", "maxQ2", "maximal Q2 value [GeV^2]", false, 80., "double", cmd);
    TCLAP::ValueArg<double> arg_minW("", "minW", "minimum W value [GeV]", false, 40., "double", cmd);
    TCLAP::ValueArg<double> arg_maxW("", "maxW", "maximal W value [GeV]", false, 180., "double", cmd);
    TCLAP::ValueArg<double> arg_minZ("", "minZ", "minimum z value", false, 0.2, "double", cmd);
    TCLAP::ValueArg<double> arg_maxZ("", "maxZ", "maximal z value", false, 0.9, "double", cmd);
    TCLAP::SwitchArg arg_gauge("g", "gauge", "test gauge invariance", cmd, false);
    TCLAP::ValueArg<string> arg_pdfName("p", "pdf", "name of PDF set", false, "cteq6l1", "string", cmd);
    TCLAP::ValueArg<int> arg_bins("b", "bins", "number histigram bins", false, 30, "int", cmd);
    TCLAP::ValueArg<double> arg_xi("x", "xi", "scale parameter", false, 1, "double", cmd);
    TCLAP::SwitchArg arg_ldme("", "ldme", "whether multiply final answers by LDMEs", cmd, false);
    TCLAP::ValueArg<string> outprefix_arg("o", "out", "prefix for outupu files", false, "", "string", cmd);
    TCLAP::SwitchArg arg_saveHist("H", "saveHist", "whether to save histograms", cmd, false);
    TCLAP::ValueArg<int> arg_seed("", "seed", "random seed", false, 0, "int", cmd);
    TCLAP::SwitchArg arg_check("c", "check", "whether to check kinematics", cmd, false);

    cmd.parse(argc, argv);
    nEv = (int) arg_nEv.getValue();
    ecm = (dbl_type) arg_ecm.getValue();
    minQ2 = (dbl_type) arg_minQ2.getValue();
    maxQ2 = (dbl_type) arg_maxQ2.getValue();
    minW2=pow((dbl_type)arg_minW.getValue(),2);
    maxW2=pow((dbl_type)arg_maxW.getValue(),2);
    zMin=arg_minZ.getValue();
    zMax=arg_maxZ.getValue();
    gauge=arg_gauge.getValue();
    pdf_set_name=arg_pdfName.getValue();
    nBins=arg_bins.getValue();
    xi=arg_xi.getValue();
    ldme=arg_ldme.getValue();
    if(!ldme) {
        O3S11=1;
        O3S18=1;
        O3P08 = 1;
        O1S08 = 1;
    };
    out_prefix=outprefix_arg.getValue();
    saveHist=arg_saveHist.getValue();
    seed=arg_seed.getValue();
    check = arg_check.getValue();
}

const int nChannels=4;
string channel_name[nChannels];
TH1F *hPT2[nChannels], *hQ2[nChannels], *hY[nChannels], *hZ[nChannels], *hW[nChannels];
dbl_type MATR2[nChannels], Sigma[nChannels];
//std::function<dbl_type (void)> matr2_func[nChannels];

const int nVars=8;
string tuple_vars="Q2:Y:pTpsi:W2:z:x:wt:pdf";
Float_t tuple_vals[nVars+nChannels];

void init_channels(void) {
    channel_name[0]="3S1_cs"; //matr2_func[0]=getMatr2_3S1_cs();
    channel_name[1]="3S1_co"; //matr2_func[1]=getMatr2_3S1_co();
    channel_name[2]="1S0_co"; //matr2_func[2]=getMatr2_1S0_co();
    channel_name[3]="3P0_co"; //matr2_func[3]=getMatr2_3P0_co();
    for(int iChannel=0; iChannel<nChannels; ++iChannel) {
        tuple_vars = tuple_vars+(":matr2_"+channel_name[iChannel]);
    }
    cout<<" tuple_vars="<<tuple_vars<<endl;
    if(saveHist) {
        for(int iChannel=0; iChannel<nChannels; ++iChannel) {
            hPT2[iChannel]=new TH1F( ("hPT2_"+channel_name[iChannel]).c_str(),"hPT2",nBins,0,30);
            hPT2[iChannel]->Sumw2();
            hQ2[iChannel]=new TH1F( ("hQ2_"+channel_name[iChannel]).c_str(),"hQ2",nBins,minQ2,maxQ2);
            hQ2[iChannel]->Sumw2();
            hZ[iChannel]=new TH1F( ("hZ_"+channel_name[iChannel]).c_str(),"hZ",nBins,0,1);
            hZ[iChannel]->Sumw2();
            hY[iChannel]=new TH1F( ("hY_"+channel_name[iChannel]).c_str(),"hY",nBins,0,1);
            hY[iChannel]->Sumw2();
            hW[iChannel]=new TH1F( ("hW_"+channel_name[iChannel]).c_str(),"hW",nBins,sqrt(minW2),sqrt(maxW2));
            hW[iChannel]->Sumw2();
            Sigma[iChannel]=0;
        };
    };
}

void fill_hst() {
//    for(int iChannel=0; iChannel<nChannels; ++iChannel) {
//        hPT2[iChannel]->Fill(pT_squared(pPsi), dSigma[iChannel]);
//        hQ2[iChannel]->Fill(Q2, dSigma[iChannel]);
//        hZ[iChannel]->Fill(z, dSigma[iChannel]);
//        hY[iChannel]->Fill(Y, dSigma[iChannel]);
//        hW[iChannel]->Fill(sqrt(W2), dSigma[iChannel]);        
//    }
}



void save_hst() {
//    for(int iChannel=0; iChannel<nChannels; ++iChannel) {
//        hPT2[iChannel]->Write(); write_histogram_to_file(*hPT2[iChannel],(out_prefix+"hPT2_"+channel_name[iChannel]+".hst").c_str());
//        hQ2[iChannel]->Write(); write_histogram_to_file(*hQ2[iChannel],(out_prefix+"hQ2_"+channel_name[iChannel]+".hst").c_str());
//        hZ[iChannel]->Write(); write_histogram_to_file(*hZ[iChannel],(out_prefix+"hZ_"+channel_name[iChannel]+".hst").c_str());
//        hY[iChannel]->Write(); write_histogram_to_file(*hY[iChannel],(out_prefix+"hY_"+channel_name[iChannel]+".hst").c_str());
//        hW[iChannel]->Write(); write_histogram_to_file(*hW[iChannel],(out_prefix+"hW_"+channel_name[iChannel]+".hst").c_str());
//    };
}

void print_sigma(bool _endl) {
    for(int iChannel=0; iChannel<nChannels; ++iChannel) {
        cout<<"sigma["<<channel_name[iChannel]<<"]="<<Sigma[iChannel]<<" ";
        if(_endl) cout<<endl;
    };
}

int main(int argc, char **argv) {
    init_from_command_line(argc, argv);
    cout << " pdfse=" << pdf_set_name << endl;
    cout << " xi=" << xi << endl;
    lhapdf_pdf = LHAPDF::mkPDF(pdf_set_name, 0);
    Random *random = new Random(seed);

    // experimental conditions from Kinehl
    dbl_type S = pow(ecm, 2);
    RamboEP *ramEP = new RamboEP(ecm, random, Mcc, 0);
    ramEP->set_minmaxQ2(minQ2, maxQ2);
    dbl_type minX = (minQ2 + minQ2) / S / 2, maxX = (2 * maxW2 + maxQ2) / S;
    if (maxX > 1) maxX = 1;
    minX = 0;
    maxX = 1;
    cout << " minX=" << minX << " maxX=" << maxX << endl;
    ramEP->set_minmaxY(minX, maxX);

    cout << " gauge: " << gauge << endl;
    cout << " ldme "<<ldme<<endl;
    cout << "out_prefix="<<out_prefix<<endl;
    cout << " saveHist="<<saveHist<<endl;
    cout << " seed=" << seed<<endl;
    TFile file((out_prefix+"_out.root").c_str(), "RECREATE");
    
    init_channels();
    TNtuple tup("tup", "tup", tuple_vars.c_str());
    TNtuple stats_tuple("stats","stats","nEv:nPassed:seed");


    dbl_type Q2_scale;
    int nPassed = 0, nNegative = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        if (iEv % (nEv / 100) == 0 && iEv > 0) {
            cout << "---- Event " << iEv << " ("<< (int) (100. * iEv / nEv) << " %)";
            if(ldme) cout<<" --- sigmaCS=" << Sigma[0] * nEv / iEv << " pb" << endl;
            else cout<<" --- sigmaCS/LDME=" << Sigma[0] * nEv / iEv << " pb*GeV^(-3)" << endl;
        };
        x = random->rand(minX, maxX);

        
        
        if (!kinematics(ramEP)) continue;
        dbl_type wt = (maxX-minX)*ramEP->wt;
        Q2 = -mass2(k1);
        set_v4(P, ramEP->P);
        Y = sp(k1, P) / sp(k, P);
        z = sp(pPsi, P) / sp(k1, P);
        W2 = sum_mass2(k1, P);


        if (W2 < minW2 || W2 > maxW2) continue;
        if (z < zMin || z > zMax) continue;
        if (Q2 < minQ2 || Q2 > maxQ2) continue;

        
        Q2_scale = pow(xi * pT(pPsi), 2);
        Q2_scale = xi * xi * (Q2 + Mcc * Mcc);
        alphas = lhapdf_pdf->alphasQ2((double) Q2_scale);
        dbl_type pdf = lhapdf_pdf->xfxQ2(0, (double) x, (double) Q2_scale) / x;

        dbl_type matr2[nChannels];
        matr2[0]=getMatr2_3S1_cs();
        matr2[1]=getMatr2_3S1_co();
        matr2[2]=getMatr2_1S0_co();
        matr2[3]=getMatr2_3P0_co();
        
//        dbl_type norm_sigma=1;
//        norm_sigma *= pdf*wt;
        wt *= 1. / 4 / 8; // polarizations and initial gluon spin
        wt *= 1. / (2 * x * pow(ecm, 2)); // 4 I
        wt *= 1. / nEv;
        wt *= picob;

        tuple_vals[0]=Q2;
        tuple_vals[1]=Y;
        tuple_vals[2]=pT(pPsi);
        tuple_vals[3]=W2;
        tuple_vals[4]=z;
        tuple_vals[5]=x;
        tuple_vals[6]=wt;
        tuple_vals[7]=pdf;
        for(int iChannel=0; iChannel<nChannels; ++iChannel) {
            MATR2[iChannel]=matr2[iChannel];
            tuple_vals[nVars+iChannel]=MATR2[iChannel]*nEv;
            Sigma[iChannel] += wt*pdf*MATR2[iChannel];
        };

        tup.Fill(tuple_vals);
        if(saveHist) fill_hst();
        ++nPassed;
    };
    stats_tuple.Fill(nEv,nPassed,seed);
    tup.Write(); stats_tuple.Write();
    if(saveHist) save_hst();
    file.Save();

    print_sigma(true);
    cout << nPassed << " (" << (int) (100. * nPassed / nEv) << "%) events passed" << endl;
    cout << nNegative << " (" << (int) (100. * nNegative / nEv) << "%) events with negative matr2" << endl;
    cout << ramEP->nFault << " faults in ramEP" << endl;

}