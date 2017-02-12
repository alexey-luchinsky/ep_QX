#include "kinematics/RamboEP.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>

#include <tclap/CmdLine.h>


dbl_type Mcc = 3.1, mc = Mcc / 2, Opsi = 0.270;
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
dbl_type minQ2, maxQ2;
dbl_type minW2, maxW2;
dbl_type zMin, zMax;
bool gauge;
string pdf_set_name;
int nBins;
dbl_type xi;





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
}

int main(int argc, char **argv) {
    init_from_command_line(argc, argv);
    cout << " pdfse=" << pdf_set_name << endl;
    cout << " xi=" << xi << endl;
    lhapdf_pdf = LHAPDF::mkPDF(pdf_set_name, 0);
    Random *random = new Random();

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

    TFile file("out_pol.root", "RECREATE");
    TNtuple tup("tup", "tup", "Q2:Y:pTpsi:W2:dsigma");
    TH1F hPT2_cs("hPT2_cs", "hPT2_cs", nBins, 0, 30);
    hPT2_cs.Sumw2();
    TH1F hQ2_cs("hQ2_cs", "hQ2_cs", nBins, minQ2, maxQ2);
    hQ2_cs.Sumw2();
    TH1F hZ_cs("hZ_cs", "hZ_cs", nBins, 0, 1);
    hZ_cs.Sumw2();
    TH1F hY_cs("hY_cs", "hY_cs", nBins, 0, 1);
    hY_cs.Sumw2();
    TH1F hW_cs("hW_cs", "hW_cs", nBins, sqrt(minW2), sqrt(maxW2));
    hW_cs.Sumw2();

    TH1F hPT2_co("hPT2_co", "hPT2_co", nBins, 0, 30);
    hPT2_co.Sumw2();
    TH1F hQ2_co("hQ2_co", "hQ2_co", nBins, minQ2, maxQ2);
    hQ2_co.Sumw2();
    TH1F hZ_co("hZ_co", "hZ_co", nBins, 0, 1);
    hZ_co.Sumw2();
    TH1F hY_co("hY_co", "hY_co", nBins, 0, 1);
    hY_co.Sumw2();
    TH1F hW_co("hW_co", "hW_co", nBins, sqrt(minW2), sqrt(maxW2));
    hW_co.Sumw2();

    TH1F hPT2_co2("hPT2_co2", "hPT2_co2", nBins, 0, 30);
    hPT2_co2.Sumw2();
    TH1F hQ2_co2("hQ2_co2", "hQ2_co2", nBins, minQ2, maxQ2);
    hQ2_co2.Sumw2();
    TH1F hZ_co2("hZ_co2", "hZ_co2", nBins, 0, 1);
    hZ_co2.Sumw2();
    TH1F hY_co2("hY_co2", "hY_co2", nBins, 0, 1);
    hY_co2.Sumw2();
    TH1F hW_co2("hW_co2", "hW_co2", nBins, sqrt(minW2), sqrt(maxW2));
    hW_co2.Sumw2();

    TH1F hPT2_co3("hPT2_co3", "hPT2_co3", nBins, 0, 30);
    hPT2_co3.Sumw2();
    TH1F hQ2_co3("hQ2_co3", "hQ2_co3", nBins, minQ2, maxQ2);
    hQ2_co3.Sumw2();
    TH1F hZ_co3("hZ_co3", "hZ_co3", nBins, 0, 1);
    hZ_co3.Sumw2();
    TH1F hY_co3("hY_co3", "hY_co3", nBins, 0, 1);
    hY_co3.Sumw2();
    TH1F hW_co3("hW_co3", "hW_co3", nBins, sqrt(minW2), sqrt(maxW2));
    hW_co3.Sumw2();


    dbl_type Q2_scale;
    int nPassed = 0, nNegative = 0;
    dbl_type sum = 0, dsigmaCS, sigmaCS = 0, dsigmaCO, sigmaCO = 0, dsigmaCO2, sigmaCO2 = 0, dsigmaCO3, sigmaCO3 = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        if (iEv % (nEv / 100) == 0 && iEv > 0) cout << "---- Event " << iEv << " ("
                << (int) (100. * iEv / nEv) << " %) --- sigmaCS=" << sigmaCS * nEv / iEv << " pb"
                " sigmaCO=" << sigmaCO * nEv / iEv << " pb" <<
                " sigmaCO2=" << sigmaCO2 * nEv / iEv << " pb" <<
                " sigmaCO3=" << sigmaCO3 * nEv / iEv << " pb" << endl;
        x = random->rand(minX, maxX);

        if (!kinematics(ramEP)) continue;
        ramEP->wt *= maxX - minX;
        dbl_type Q2 = -mass2(k1);
        set_v4(P, ramEP->P);
        dbl_type Y = sp(k1, P) / sp(k, P);
        dbl_type z = sp(pPsi, P) / sp(k1, P);
        dbl_type W2 = sum_mass2(k1, P);


        if (W2 < minW2 || W2 > maxW2) continue;
        if (z < zMin || z > zMax) continue;
        if (Q2 < minQ2 || Q2 > maxQ2) continue;

        Q2_scale = pow(xi * pT(pPsi), 2);
        Q2_scale = xi * xi * (Q2 + Mcc * Mcc);
        alphas = lhapdf_pdf->alphasQ2((double) Q2_scale);
        dbl_type pdf = lhapdf_pdf->xfxQ2(0, (double) x, (double) Q2_scale) / x;

        dbl_type wt = ramEP->wt;
        dbl_type matr2CS = getMatr2_3S1_cs();
        dbl_type matr2CO = getMatr2_3S1_co();
        dbl_type matr2CO2 = getMatr2_1S0_co();
        dbl_type matr2CO3 = getMatr2_3P0_co();
        if (!(matr2CS > 0 && matr2CO > 0 && matr2CO2 > 0)) {
            nNegative++;
            continue;
        };
        
        dbl_type norm_sigma=1;
        norm_sigma *= pdf*wt;
        norm_sigma *= 1. / 4 / 8; // polarizations and initial gluon spin
        norm_sigma *= 1. / (2 * x * pow(ecm, 2)); // 4 I
        norm_sigma *= 1. / nEv;
        norm_sigma *= picob;
        
        dsigmaCS = matr2CS * norm_sigma;   sigmaCS += dsigmaCS;
        dsigmaCO = matr2CO * norm_sigma;   sigmaCO += dsigmaCO;
        dsigmaCO2 = matr2CO2 * norm_sigma; sigmaCO2 += dsigmaCO2;
        dsigmaCO3 = matr2CO3 * norm_sigma; sigmaCO3 += dsigmaCO3;

        tup.Fill(ramEP->Q2, ramEP->Y, pT(pPsi), W2, dsigmaCS);
        hPT2_cs.Fill(pT_squared(pPsi), dsigmaCS);
        hQ2_cs.Fill(ramEP->Q2, dsigmaCS);
        hZ_cs.Fill(z, dsigmaCS);
        hY_cs.Fill(ramEP->Y, dsigmaCS);
        hW_cs.Fill(sqrt(W2), dsigmaCS);

        hPT2_co.Fill(pT_squared(pPsi), dsigmaCO);
        hQ2_co.Fill(ramEP->Q2, dsigmaCO);
        hZ_co.Fill(z, dsigmaCO);
        hY_co.Fill(ramEP->Y, dsigmaCO);
        hW_co.Fill(sqrt(W2), dsigmaCO);

        hPT2_co2.Fill(pT_squared(pPsi), dsigmaCO2);
        hQ2_co2.Fill(ramEP->Q2, dsigmaCO2);
        hZ_co2.Fill(z, dsigmaCO2);
        hY_co2.Fill(ramEP->Y, dsigmaCO2);
        hW_co2.Fill(sqrt(W2), dsigmaCO2);

        hPT2_co3.Fill(pT_squared(pPsi), dsigmaCO3);
        hQ2_co3.Fill(ramEP->Q2, dsigmaCO3);
        hZ_co3.Fill(z, dsigmaCO3);
        hY_co3.Fill(ramEP->Y, dsigmaCO3);
        hW_co3.Fill(sqrt(W2), dsigmaCO3);

        ++nPassed;
    };

    tup.Write();
    hPT2_cs.Write();
    file.Save();
    write_histogram_to_file(hPT2_cs, "hPT2_cs.hst");
    write_histogram_to_file(hQ2_cs, "hQ2_cs.hst");
    write_histogram_to_file(hZ_cs, "hZ_cs.hst");
    write_histogram_to_file(hY_cs, "hY_cs.hst");
    write_histogram_to_file(hW_cs, "hW_cs.hst");

    write_histogram_to_file(hPT2_co, "hPT2_co.hst");
    write_histogram_to_file(hQ2_co, "hQ2_co.hst");
    write_histogram_to_file(hZ_co, "hZ_co.hst");
    write_histogram_to_file(hY_co, "hY_co.hst");
    write_histogram_to_file(hW_co, "hW_co.hst");

    write_histogram_to_file(hPT2_co2, "hPT2_co2.hst");
    write_histogram_to_file(hQ2_co2, "hQ2_co2.hst");
    write_histogram_to_file(hZ_co2, "hZ_co2.hst");
    write_histogram_to_file(hY_co2, "hY_co2.hst");
    write_histogram_to_file(hW_co2, "hW_co2.hst");

    write_histogram_to_file(hPT2_co3, "hPT2_co3.hst");
    write_histogram_to_file(hQ2_co3, "hQ2_co3.hst");
    write_histogram_to_file(hZ_co3, "hZ_co3.hst");
    write_histogram_to_file(hY_co3, "hY_co3.hst");
    write_histogram_to_file(hW_co3, "hW_co3.hst");


    cout << " sigmaCS=" << sigmaCS << " pb" << endl;
    cout << " sigmaCO=" << sigmaCO << " pb" << endl;
    cout << " sigmaCO2=" << sigmaCO2 << " pb" << endl;
    cout << nPassed << " (" << (int) (100. * nPassed / nEv) << "%) events passed" << endl;
    cout << nNegative << " (" << (int) (100. * nNegative / nEv) << "%) events with negative matr2" << endl;
    cout << ramEP->nFault << " faults in ramEP" << endl;

}