#include "kinematics/RamboEP.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>

dbl_type Mcc=3.1, mc=Mcc/2, Opsi=0.270;
dbl_type nanob=0.389e6, picob=1e3*nanob; // conversion to barn
dbl_type PI=acos(-1), alpha=1./137, alphas=0.3;
dbl_type ecm=10, x;
bool gauge=true;

#include "LHAPDF/LHAPDF.h"
LHAPDF::PDF *lhapdf_pdf;

void write_histogram_to_file(TH1F &histogram, string file_name);
bool kinematics(RamboEP *ramEP);
dbl_type getMatr2_3S1_cs(dbl_type(&k)[4], dbl_type(&kp)[4], dbl_type(&k1)[4], dbl_type(&k2)[4], 
                dbl_type(&k3)[4], dbl_type(&pPsi)[4], int iEv);
dbl_type getMatr2_3S1_co(dbl_type(&k)[4], dbl_type(&kp)[4], dbl_type(&k1)[4], dbl_type(&k2)[4], 
                dbl_type(&k3)[4], dbl_type(&pPsi)[4], int iEv);
dbl_type getMatr2_1S0_co(dbl_type(&k)[4], dbl_type(&kp)[4], dbl_type(&k1)[4], dbl_type(&k2)[4], 
                dbl_type(&k3)[4], dbl_type(&pPsi)[4], int iEv);
dbl_type getMatr2_3P0_co(dbl_type(&k)[4], dbl_type(&kp)[4], dbl_type(&k1)[4], dbl_type(&k2)[4], 
                dbl_type(&k3)[4], dbl_type(&pPsi)[4], int iEv);

dbl_type P[4], k[4], kp[4], k1[4], k2[4], k3[4], pPsi[4];


int main(void) {
    string pdf_set_name="cteq6l1"; cout<<" pdfse="<<pdf_set_name<<endl;
    dbl_type xi=1; cout<<" xi="<<xi<<endl;
    lhapdf_pdf=LHAPDF::mkPDF(pdf_set_name,0);
    Random *random=new Random();

    // experimental conditions from Kinehl
    ecm=300;
    cout<<" ecm="<<ecm<<endl;
    dbl_type S=pow(ecm,2);
    RamboEP *ramEP=new RamboEP(ecm, random, Mcc,0);
    dbl_type minQ2=4, maxQ2=80;
    ramEP->set_minmaxQ2(minQ2, maxQ2);
    dbl_type minW2=pow(40,2), maxW2=pow(180,2);
    dbl_type minX=(minQ2+minQ2)/S/2, maxX=(2*maxW2+maxQ2)/S;
    if(maxX>1) maxX=1;
    minX=0; maxX=1;
    cout<<" minX="<<minX<<" maxX="<<maxX<<endl;
    ramEP->set_minmaxY(minX,maxX);
    dbl_type zMin=0.2, zMax=0.9; 

    gauge=true;
    cout<<" gauge: "<<gauge<<endl;
    
    int nBins=30;
    TFile file("out_pol.root","RECREATE");
    TNtuple tup("tup","tup","Q2:Y:pTpsi:W2:dsigma");
    TH1F hPT2_cs("hPT2_cs","hPT2_cs",nBins,0,30); hPT2_cs.Sumw2();
    TH1F hQ2_cs("hQ2_cs","hQ2_cs",nBins,minQ2,maxQ2); hQ2_cs.Sumw2();
    TH1F hZ_cs("hZ_cs","hZ_cs",nBins,0,1); hZ_cs.Sumw2();
    TH1F hY_cs("hY_cs","hY_cs",nBins,0,1); hY_cs.Sumw2();
    TH1F hW_cs("hW_cs","hW_cs",nBins,sqrt(minW2),sqrt(maxW2)); hW_cs.Sumw2();
    
    TH1F hPT2_co("hPT2_co","hPT2_co",nBins,0,30); hPT2_co.Sumw2();
    TH1F hQ2_co("hQ2_co","hQ2_co",nBins,minQ2,maxQ2); hQ2_co.Sumw2();
    TH1F hZ_co("hZ_co","hZ_co",nBins,0,1); hZ_co.Sumw2();
    TH1F hY_co("hY_co","hY_co",nBins,0,1); hY_co.Sumw2();
    TH1F hW_co("hW_co","hW_co",nBins,sqrt(minW2),sqrt(maxW2)); hW_co.Sumw2();
    
    TH1F hPT2_co2("hPT2_co2","hPT2_co2",nBins,0,30); hPT2_co2.Sumw2();
    TH1F hQ2_co2("hQ2_co2","hQ2_co2",nBins,minQ2,maxQ2); hQ2_co2.Sumw2();
    TH1F hZ_co2("hZ_co2","hZ_co2",nBins,0,1); hZ_co2.Sumw2();
    TH1F hY_co2("hY_co2","hY_co2",nBins,0,1); hY_co2.Sumw2();
    TH1F hW_co2("hW_co2","hW_co2",nBins,sqrt(minW2),sqrt(maxW2)); hW_co2.Sumw2();
    
    TH1F hPT2_co3("hPT2_co3","hPT2_co3",nBins,0,30); hPT2_co3.Sumw2();
    TH1F hQ2_co3("hQ2_co3","hQ2_co3",nBins,minQ2,maxQ2); hQ2_co3.Sumw2();
    TH1F hZ_co3("hZ_co3","hZ_co3",nBins,0,1); hZ_co3.Sumw2();
    TH1F hY_co3("hY_co3","hY_co3",nBins,0,1); hY_co3.Sumw2();
    TH1F hW_co3("hW_co3","hW_co3",nBins,sqrt(minW2),sqrt(maxW2)); hW_co3.Sumw2();


    dbl_type Q2_scale;
    int nEv=1e7, nPassed=0, nNegative=0;
    dbl_type sum=0, dsigmaCS, sigmaCS=0, dsigmaCO, sigmaCO=0, dsigmaCO2, sigmaCO2=0, dsigmaCO3, sigmaCO3=0;
    for(int iEv=0; iEv<nEv; ++iEv) {
        if( iEv % (nEv/100) == 0 && iEv>0) cout<<"---- Event "<<iEv<<" ("
                <<(int)(100.*iEv/nEv)<<" %) --- sigmaCS="<<sigmaCS*nEv/iEv<<" pb"
                " sigmaCO="<<sigmaCO*nEv/iEv<<" pb"<<
                " sigmaCO2="<<sigmaCO2*nEv/iEv<<" pb"<<
                " sigmaCO3="<<sigmaCO3*nEv/iEv<<" pb"<<endl;
        x=random->rand(minX,maxX);

        if(!kinematics(ramEP)) continue;
        ramEP->wt *= maxX-minX;
        dbl_type Q2=-mass2(k1);
        set_v4(P,ramEP->P);
        dbl_type Y=sp(k1,P)/sp(k,P);
        dbl_type z=sp(pPsi,P)/sp(k1,P);
        dbl_type W2=sum_mass2(k1,P);
        
//        dbl_type q[4];
//        subtract(k,kp,q);
//        assert(is_zero(spp(k,kp,k1,q)));       

        if(W2<minW2||W2>maxW2) continue;
        if(z<zMin || z>zMax) continue;
        if(Q2<minQ2||Q2>maxQ2) continue;
            
        Q2_scale=pow(xi*pT(pPsi),2);
        Q2_scale = xi*xi*(Q2+Mcc*Mcc);
        alphas=lhapdf_pdf->alphasQ2((double) Q2_scale);
        dbl_type pdf = lhapdf_pdf->xfxQ2(0, (double)x, (double)Q2_scale)/x;

        dbl_type wt=ramEP->wt;
        dbl_type matr2CS=getMatr2_3S1_cs(k,kp,k1,k2,k3,pPsi,nPassed);
        dbl_type matr2CO=getMatr2_3S1_co(k,kp,k1,k2,k3,pPsi,nPassed);
        dbl_type matr2CO2=getMatr2_1S0_co(k,kp,k1,k2,k3,pPsi,nPassed);
        dbl_type matr2CO3=getMatr2_3P0_co(k,kp,k1,k2,k3,pPsi,nPassed);
        if(!(matr2CS>0 && matr2CO>0 && matr2CO2>0)) {
            nNegative++;
            continue;
        };      

        dsigmaCS=1;
        dsigmaCS *= matr2CS*pdf*wt;
        dsigmaCS *= 1./4/8; // polarizations and initial gluon spin
        dsigmaCS *= 1./(2*x*pow(ecm,2)); // 4 I
        dsigmaCS *= 1./nEv;
        dsigmaCS *= picob;
        sigmaCS += dsigmaCS;

        dsigmaCO=1;
        dsigmaCO *= matr2CO*pdf*wt;
        dsigmaCO *= 1./4/8; // polarizations and initial gluon spin
        dsigmaCO *= 1./(2*x*pow(ecm,2)); // 4 I
        dsigmaCO *= 1./nEv;
        dsigmaCO *= picob;
        sigmaCO += dsigmaCO;

        dsigmaCO2=1;
        dsigmaCO2 *= matr2CO2*pdf*wt;
        dsigmaCO2 *= 1./4/8; // polarizations and initial gluon spin
        dsigmaCO2 *= 1./(2*x*pow(ecm,2)); // 4 I
        dsigmaCO2 *= 1./nEv;
        dsigmaCO2 *= picob;
        sigmaCO2 += dsigmaCO2;
        
        dsigmaCO3=1;
        dsigmaCO3 *= matr2CO3*pdf*wt;
        dsigmaCO3 *= 1./4/8; // polarizations and initial gluon spin
        dsigmaCO3 *= 1./(2*x*pow(ecm,2)); // 4 I
        dsigmaCO3 *= 1./nEv;
        dsigmaCO3 *= picob;
        sigmaCO3 += dsigmaCO3;

        tup.Fill(ramEP->Q2,ramEP->Y,pT(pPsi), W2, dsigmaCS);
        hPT2_cs.Fill(pT_squared(pPsi),dsigmaCS);
        hQ2_cs.Fill(ramEP->Q2,dsigmaCS);
        hZ_cs.Fill(z,dsigmaCS);
        hY_cs.Fill(ramEP->Y,dsigmaCS);
        hW_cs.Fill(sqrt(W2),dsigmaCS);
        
        hPT2_co.Fill(pT_squared(pPsi),dsigmaCO);
        hQ2_co.Fill(ramEP->Q2,dsigmaCO);
        hZ_co.Fill(z,dsigmaCO);
        hY_co.Fill(ramEP->Y,dsigmaCO);
        hW_co.Fill(sqrt(W2),dsigmaCO);

        hPT2_co2.Fill(pT_squared(pPsi),dsigmaCO2);
        hQ2_co2.Fill(ramEP->Q2,dsigmaCO2);
        hZ_co2.Fill(z,dsigmaCO2);
        hY_co2.Fill(ramEP->Y,dsigmaCO2);
        hW_co2.Fill(sqrt(W2),dsigmaCO2);

        hPT2_co3.Fill(pT_squared(pPsi),dsigmaCO3);
        hQ2_co3.Fill(ramEP->Q2,dsigmaCO3);
        hZ_co3.Fill(z,dsigmaCO3);
        hY_co3.Fill(ramEP->Y,dsigmaCO3);
        hW_co3.Fill(sqrt(W2),dsigmaCO3);

        ++nPassed;
    };
    
    tup.Write(); hPT2_cs.Write();
    file.Save();
    write_histogram_to_file(hPT2_cs,"hPT2_cs.hst");
    write_histogram_to_file(hQ2_cs,"hQ2_cs.hst");
    write_histogram_to_file(hZ_cs,"hZ_cs.hst");
    write_histogram_to_file(hY_cs,"hY_cs.hst");
    write_histogram_to_file(hW_cs,"hW_cs.hst");

    write_histogram_to_file(hPT2_co,"hPT2_co.hst");
    write_histogram_to_file(hQ2_co,"hQ2_co.hst");
    write_histogram_to_file(hZ_co,"hZ_co.hst");
    write_histogram_to_file(hY_co,"hY_co.hst");
    write_histogram_to_file(hW_co,"hW_co.hst");

    write_histogram_to_file(hPT2_co2,"hPT2_co2.hst");
    write_histogram_to_file(hQ2_co2,"hQ2_co2.hst");
    write_histogram_to_file(hZ_co2,"hZ_co2.hst");
    write_histogram_to_file(hY_co2,"hY_co2.hst");
    write_histogram_to_file(hW_co2,"hW_co2.hst");

    write_histogram_to_file(hPT2_co3,"hPT2_co3.hst");
    write_histogram_to_file(hQ2_co3,"hQ2_co3.hst");
    write_histogram_to_file(hZ_co3,"hZ_co3.hst");
    write_histogram_to_file(hY_co3,"hY_co3.hst");
    write_histogram_to_file(hW_co3,"hW_co3.hst");


    cout<<" sigmaCS="<<sigmaCS<<" pb"<<endl;
    cout<<" sigmaCO="<<sigmaCO<<" pb"<<endl;
    cout<<" sigmaCO2="<<sigmaCO2<<" pb"<<endl;
    cout<<nPassed<<" ("<<(int)(100.*nPassed/nEv)<<"%) events passed"<<endl;
    cout<<nNegative<<" ("<<(int)(100.*nNegative/nEv)<<"%) events with negative matr2"<<endl;
    cout<<ramEP->nFault<<" faults in ramEP"<<endl;

 }