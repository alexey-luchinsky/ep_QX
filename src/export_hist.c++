#include "algebra.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TChain.h"
using namespace std;

void write_histogram_to_file(TH1F &histogram, string file_name) {
    const char *__file_name__ = file_name.c_str();
    remove(__file_name__);
    ofstream file;

    file.open(__file_name__);
    for (int i = 1; i <= histogram.GetNbinsX(); i++)
        file << setiosflags(ios::scientific) << histogram.GetBinCenter(i) <<
        " " << setiosflags(ios::scientific) << histogram.GetBinContent(i) / histogram.GetBinWidth(i) <<
        " " << setiosflags(ios::scientific) << histogram.GetBinError(i) / histogram.GetBinWidth(i) << endl;

    file.close();
}
const int nVars=5;
const int nChannels=4;
string channel_name[nChannels];
TH1F *hPT2[nChannels], *hQ2[nChannels], *hY[nChannels], *hZ[nChannels], *hW[nChannels];
dbl_type dSigma[nChannels], Sigma[nChannels];


string tuple_vars[nVars+nChannels];
Float_t tuple_vals[nVars+nChannels];

TChain tup_chain("tup"), stat_chain("stats");


void init_channels(void) {
    tuple_vars[0]="Q2"; 
    tuple_vars[1]="Y";
    tuple_vars[2]="pTpsi";
    tuple_vars[3]="W2";
    tuple_vars[4]="z";
    channel_name[0]="3S1_cs";
    channel_name[1]="3S1_co";
    channel_name[2]="1S0_co";
    channel_name[3]="3P0_co";
    for(int iChannel=0; iChannel<nChannels; ++iChannel)
        tuple_vars[nVars+iChannel]="s"+channel_name[iChannel];
    for(int i=0; i<nVars+nChannels; ++i) {
        tup_chain.SetBranchAddress(tuple_vars[i].c_str(), &(tuple_vals[i]));
    };
}

int nBins=30;
string out_prefix="";
dbl_type Q2, minQ2, maxQ2, Y, pTpsi, minPT2, maxPT2, W2, minW2, maxW2, z;
int nEv=0;

void init_hists() {
    for(int iChannel=0; iChannel<nChannels; ++iChannel) {
        minPT2=0;
        maxPT2=30;
        hPT2[iChannel]=new TH1F( ("hPT2_"+channel_name[iChannel]).c_str(),"hPT2",nBins,minPT2,maxPT2);
        hPT2[iChannel]->Sumw2();
        //
        minQ2=tup_chain.GetMinimum("Q2"); maxQ2=tup_chain.GetMaximum("Q2");
        hQ2[iChannel]=new TH1F( ("hQ2_"+channel_name[iChannel]).c_str(),"hQ2",nBins,minQ2, maxQ2);
        hQ2[iChannel]->Sumw2();
        //
        hZ[iChannel]=new TH1F( ("hZ_"+channel_name[iChannel]).c_str(),"hZ",nBins,0,1);
        hZ[iChannel]->Sumw2();
        hY[iChannel]=new TH1F( ("hY_"+channel_name[iChannel]).c_str(),"hY",nBins,0,1);
        hY[iChannel]->Sumw2();
        minW2=tup_chain.GetMinimum("W2"); maxW2=tup_chain.GetMaximum("W2");
        hW[iChannel]=new TH1F( ("hW_"+channel_name[iChannel]).c_str(),"hW",nBins,sqrt(minW2),sqrt(maxW2));
        hW[iChannel]->Sumw2();
        Sigma[iChannel]=0;
    };
    
}

void save_hst() {
    for(int iChannel=0; iChannel<nChannels; ++iChannel) {
        write_histogram_to_file(*hPT2[iChannel],(out_prefix+"hPT2_"+channel_name[iChannel]+".hst").c_str());
        write_histogram_to_file(*hQ2[iChannel],(out_prefix+"hQ2_"+channel_name[iChannel]+".hst").c_str());
        write_histogram_to_file(*hZ[iChannel],(out_prefix+"hZ_"+channel_name[iChannel]+".hst").c_str());
        write_histogram_to_file(*hY[iChannel],(out_prefix+"hY_"+channel_name[iChannel]+".hst").c_str());
        write_histogram_to_file(*hW[iChannel],(out_prefix+"hW_"+channel_name[iChannel]+".hst").c_str());
    };
}

void fill_hst() {
    for(int iChannel=0; iChannel<nChannels; ++iChannel)
        dSigma[iChannel]=tuple_vals[nVars+iChannel]/nEv;
    for(int iChannel=0; iChannel<nChannels; ++iChannel) {
        hPT2[iChannel]->Fill(pow(pTpsi,2), dSigma[iChannel]);
        hQ2[iChannel]->Fill(Q2, dSigma[iChannel]);
        hZ[iChannel]->Fill(z, dSigma[iChannel]);
        hY[iChannel]->Fill(Y, dSigma[iChannel]);
        hW[iChannel]->Fill(sqrt(W2), dSigma[iChannel]);        
    }
}


int main(void) {
    tup_chain.Add("*.root"); stat_chain.Add("*.root");

    // calculate total number of events
    float nev;
    stat_chain.SetBranchAddress("nEv",&nev);
    for(int i=0; i<stat_chain.GetEntries(); ++i) {
        stat_chain.GetEntry(i);
        cout<<" nev="<<nev<<endl;
        nEv += (int)nev;
    }
    
    init_channels();
    init_hists();
    cout<<"[tup_chain]="<<tup_chain.GetEntries()<<endl;
    cout<<"[stat_chain]="<<stat_chain.GetEntries()<<endl;
    for(int iEv=0; iEv<tup_chain.GetEntries(); ++iEv) {
        if (iEv % (nEv / 100) == 0 && iEv > 0) 
            cout << "---- Event " << iEv << " ("<< (int) (100. * iEv / tup_chain.GetEntries()) << " %)"<<endl;
        tup_chain.GetEntry(iEv);
        Q2=tuple_vals[0]; Y=tuple_vals[1]; pTpsi=tuple_vals[2]; W2=tuple_vals[3];
        z=tuple_vals[4];
        fill_hst();
    }
    save_hst();
    return 0;
}
