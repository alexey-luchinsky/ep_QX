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

#include "LHAPDF/LHAPDF.h"
LHAPDF::PDF *lhapdf_pdf;

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


#define spp(k,kp,v1,v2) 8*4*(sp(k,v1)*sp(kp,v2)-sp(v1,v2)*sp(k,kp)+sp(kp,v1)*sp(k,v2))


bool kinematics(RamboEP *ramEP, 
        dbl_type(&k)[4], dbl_type(&kp)[4], dbl_type(&k1)[4], dbl_type(&k2)[4], 
                dbl_type(&k3)[4], dbl_type(&p)[4]) {
    if(!ramEP->next(kp,p,k3,x)) return false;
    set_v4(k,ramEP->kIn);                // k
        assert(is_zero(mass2(k)));       
    set_v4(kp, ramEP->kOut);             // kp
        assert(is_zero(mass2(kp)));
    set_v4(k2,ramEP->Pg);                // k2
        assert(are_equal(sum_mass2(p,k3),ramEP->W2));
    subtract(k,kp,k1);                  // k1=k-kp
        assert(are_equal(mass2(k1),-ramEP->Q2)); 
        assert(are_equal(mass2(p),Mcc*Mcc)); 
        assert(is_zero(mass2(k3)));
        assert(are_equal(sum_mass2(k,k2),x*ecm*ecm));
        assert(are_equal(sum_mass2(kp,p,k3),x*ecm*ecm));
    return true;
};



dbl_type getMatr2_pol(dbl_type(&k)[4], dbl_type(&kp)[4], dbl_type(&k1)[4], dbl_type(&k2)[4], 
                dbl_type(&k3)[4], dbl_type(&pPsi)[4], int iEv) {
    const int nAmps=9;
    cmplx Amp[nAmps];

    cmplx multTable[nAmps][nAmps];
    cmplx epsG2[4], epsG3[4], epsPsi[4];
    cmplx cepsG2[4], cepsG3[4], cepsPsi[4];
    cmplx II=COMPLEX_ONE;

    cmplx imatr2=COMPLEX_ZERO;
     for(int iPolPsi=0; iPolPsi<3; ++iPolPsi) {
         for(int iPolG2=0; iPolG2<2; ++iPolG2) {
             for(int iPolG3=0; iPolG3<2; ++iPolG3) {
                 set_gluon_polarization(iPolG2, k2, epsG2, cepsG2);
//                 for(int _i=0; _i<4; ++_i) epsG2[_i]=k2[_i];
                 set_gluon_polarization(iPolG3, k3, epsG3, cepsG3);
                 set_psi_polarization(iPolPsi, pPsi, epsPsi, cepsPsi);
                 
                dbl_type Q2=-mass2(k1);
                Amp[0]=(-64*alpha*alphas*pow(PI,2)*(pow(sp(k2,k3),2)*sp(epsG2,epsG3) - sp(epsG2,pPsi)*sp(epsG3,k3)*sp(k2,k3) - sp(epsG2,pPsi)*sp(epsG3,pPsi)*sp(k2,k3) + sp(epsG2,epsG3)*sp(k2,k3)*sp(k2,pPsi) + sp(epsG2,k3)*((sp(epsG3,k3) + sp(epsG3,pPsi))*sp(k2,pPsi) - sp(epsG3,k2)*(sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))) + sp(epsG2,pPsi)*sp(epsG3,k2)*sp(k3,pPsi) - sp(epsG2,epsG3)*sp(k2,k3)*sp(k3,pPsi) - sp(epsG2,epsG3)*sp(k2,pPsi)*sp(k3,pPsi) + sp(epsG2,k2)*(sp(epsG3,k3)*sp(k2,k3) + sp(epsG3,pPsi)*sp(k2,k3) - sp(epsG3,k2)*sp(k3,pPsi)))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
                Amp[1]=(-16*alpha*alphas*pow(PI,2)*(4*sp(epsG2,k3)*sp(epsG3,epsPsi)*(sp(k2,k3) - sp(k3,pPsi)) + 4*sp(epsG2,epsPsi)*(sp(epsG3,k3)*sp(k2,k3) + sp(epsG3,pPsi)*sp(k2,k3) - sp(epsG3,k2)*sp(k3,pPsi)) + sp(epsG2,epsG3)*(-4*sp(epsPsi,k3)*(sp(k2,k3) - sp(k3,pPsi)) + sp(epsPsi,pPsi)*(-2*sp(k2,k3) + sp(k3,pPsi))))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
                Amp[2]=(-16*alpha*alphas*pow(PI,2)*(-4*pow(sp(k2,k3),2)*sp(epsG3,epsPsi) + 4*sp(epsG3,k2)*sp(epsPsi,k3)*sp(k2,k3) + 2*sp(epsG3,k2)*sp(epsPsi,pPsi)*sp(k2,k3) + sp(epsG3,k3)*(-4*sp(epsPsi,k2)*sp(k2,k3) + sp(epsPsi,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi))) + sp(epsG3,pPsi)*(-4*sp(epsPsi,k2)*sp(k2,k3) + sp(epsPsi,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi))) + 4*sp(epsG3,k2)*sp(epsPsi,k2)*sp(k3,pPsi) - 4*sp(epsG3,k2)*sp(epsPsi,k3)*sp(k3,pPsi) - 3*sp(epsG3,k2)*sp(epsPsi,pPsi)*sp(k3,pPsi) + 4*sp(epsG3,epsPsi)*sp(k2,k3)*sp(k3,pPsi))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
                Amp[3]=(16*alpha*alphas*pow(PI,2)*(4*pow(sp(k2,k3),2)*sp(epsG2,epsPsi) - 4*sp(epsG2,pPsi)*sp(epsPsi,k3)*sp(k2,k3) - 2*sp(epsG2,pPsi)*sp(epsPsi,pPsi)*sp(k2,k3) + 4*sp(epsG2,epsPsi)*sp(k2,k3)*sp(k2,pPsi) + sp(epsG2,k3)*(4*sp(epsPsi,k3)*sp(k2,pPsi) - 4*sp(epsPsi,k2)*(sp(k2,k3) + sp(k2,pPsi)) + sp(epsPsi,pPsi)*(2*sp(k2,k3) + 3*sp(k2,pPsi))) + sp(epsG2,k2)*(4*sp(epsPsi,k3)*sp(k2,k3) + sp(epsPsi,pPsi)*(2*sp(k2,k3) - sp(k3,pPsi))) + sp(epsG2,pPsi)*sp(epsPsi,pPsi)*sp(k3,pPsi))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
                Amp[4]=(-16*alpha*alphas*pow(PI,2)*(4*sp(epsG2,k2)*sp(epsG3,epsPsi)*sp(k2,k3) - 4*sp(epsG2,pPsi)*sp(epsG3,epsPsi)*sp(k2,k3) + 4*sp(epsG2,epsPsi)*sp(epsG3,k2)*sp(k2,k3) - 4*sp(epsG2,epsG3)*sp(epsPsi,k2)*sp(k2,k3) + 2*sp(epsG2,epsG3)*sp(epsPsi,pPsi)*sp(k2,k3) + 4*sp(epsG2,k3)*sp(epsG3,epsPsi)*sp(k2,pPsi) + 4*sp(epsG2,epsPsi)*sp(epsG3,k2)*sp(k2,pPsi) - 4*sp(epsG2,epsG3)*sp(epsPsi,k2)*sp(k2,pPsi) + sp(epsG2,epsG3)*sp(epsPsi,pPsi)*sp(k2,pPsi))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
                Amp[5]=(-16*alpha*alphas*II*pow(PI,2)*sp(epsG2,epsG3)*sp(epsPsi,pPsi)*(2*pow(sp(k2,k3),2) - sp(k2,k3)*sp(k3,pPsi) + sp(k2,pPsi)*sp(k3,pPsi))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,k3)*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
                Amp[6]=(-16*alpha*alphas*II*pow(PI,2)*sp(epsG2,epsG3)*sp(epsPsi,pPsi)*(2*pow(sp(k2,k3),2) + sp(k2,k3)*sp(k2,pPsi) + sp(k2,pPsi)*sp(k3,pPsi))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,k3)*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
                Amp[7]=(-16*alpha*alphas*II*pow(PI,2)*sp(epsPsi,pPsi)*(sp(epsG2,pPsi)*sp(k2,k3)*(-2*sp(k2,k3) + sp(k3,pPsi)) + sp(epsG2,k3)*(-2*pow(sp(k2,k3),2) + sp(k2,k3)*sp(k2,pPsi) - 2*sp(k2,pPsi)*sp(k3,pPsi)) + sp(epsG2,k2)*(2*pow(sp(k2,k3),2) - sp(k2,k3)*sp(k3,pPsi) + sp(k2,pPsi)*sp(k3,pPsi)))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,k3)*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
                Amp[8]=(-16*alpha*alphas*II*pow(PI,2)*sp(epsPsi,pPsi)*(sp(epsG3,pPsi)*sp(k2,k3)*(2*sp(k2,k3) + sp(k2,pPsi)) + sp(epsG3,k3)*(2*pow(sp(k2,k3),2) + sp(k2,k3)*sp(k2,pPsi) + sp(k2,pPsi)*sp(k3,pPsi)) - sp(epsG3,k2)*(2*pow(sp(k2,k3),2) + sp(k2,k3)*sp(k3,pPsi) + 2*sp(k2,pPsi)*sp(k3,pPsi)))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,k3)*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
                 

                multTable[0][0]=(160*(sp(cepsPsi,kp)*sp(epsPsi,k) + sp(cepsPsi,k)*sp(epsPsi,kp) - sp(cepsPsi,epsPsi)*sp(k,kp)))/3.;
                multTable[0][1]=(160*(sp(epsPsi,kp)*sp(k,k2) - sp(epsPsi,k2)*sp(k,kp) + sp(epsPsi,k)*sp(k2,kp)))/3.;
                multTable[0][2]=(160*(sp(cepsG2,kp)*sp(epsPsi,k) + sp(cepsG2,k)*sp(epsPsi,kp) - sp(cepsG2,epsPsi)*sp(k,kp)))/3.;
                multTable[0][3]=(160*(sp(cepsG3,kp)*sp(epsPsi,k) + sp(cepsG3,k)*sp(epsPsi,kp) - sp(cepsG3,epsPsi)*sp(k,kp)))/3.;
                multTable[0][4]=(160*(sp(epsPsi,kp)*sp(k,k3) - sp(epsPsi,k3)*sp(k,kp) + sp(epsPsi,k)*sp(k3,kp)))/3.;
                multTable[0][5]=0;
                multTable[0][6]=0;
                multTable[0][7]=0;
                multTable[0][8]=0;
                multTable[1][0]=(160*(sp(cepsPsi,kp)*sp(k,k2) - sp(cepsPsi,k2)*sp(k,kp) + sp(cepsPsi,k)*sp(k2,kp)))/3.;
                multTable[1][1]=(320*sp(k,k2)*sp(k2,kp))/3.;
                multTable[1][2]=(160*(sp(cepsG2,kp)*sp(k,k2) - sp(cepsG2,k2)*sp(k,kp) + sp(cepsG2,k)*sp(k2,kp)))/3.;
                multTable[1][3]=(160*(sp(cepsG3,kp)*sp(k,k2) - sp(cepsG3,k2)*sp(k,kp) + sp(cepsG3,k)*sp(k2,kp)))/3.;
                multTable[1][4]=(160*(-(sp(k,kp)*sp(k2,k3)) + sp(k,k3)*sp(k2,kp) + sp(k,k2)*sp(k3,kp)))/3.;
                multTable[1][5]=0;
                multTable[1][6]=0;
                multTable[1][7]=0;
                multTable[1][8]=0;
                multTable[2][0]=(160*(sp(cepsPsi,kp)*sp(epsG2,k) + sp(cepsPsi,k)*sp(epsG2,kp) - sp(cepsPsi,epsG2)*sp(k,kp)))/3.;
                multTable[2][1]=(160*(sp(epsG2,kp)*sp(k,k2) - sp(epsG2,k2)*sp(k,kp) + sp(epsG2,k)*sp(k2,kp)))/3.;
                multTable[2][2]=(160*(sp(cepsG2,kp)*sp(epsG2,k) + sp(cepsG2,k)*sp(epsG2,kp) - sp(cepsG2,epsG2)*sp(k,kp)))/3.;
                multTable[2][3]=(160*(sp(cepsG3,kp)*sp(epsG2,k) + sp(cepsG3,k)*sp(epsG2,kp) - sp(cepsG3,epsG2)*sp(k,kp)))/3.;
                multTable[2][4]=(160*(sp(epsG2,kp)*sp(k,k3) - sp(epsG2,k3)*sp(k,kp) + sp(epsG2,k)*sp(k3,kp)))/3.;
                multTable[2][5]=0;
                multTable[2][6]=0;
                multTable[2][7]=0;
                multTable[2][8]=0;
                multTable[3][0]=(160*(sp(cepsPsi,kp)*sp(epsG3,k) + sp(cepsPsi,k)*sp(epsG3,kp) - sp(cepsPsi,epsG3)*sp(k,kp)))/3.;
                multTable[3][1]=(160*(sp(epsG3,kp)*sp(k,k2) - sp(epsG3,k2)*sp(k,kp) + sp(epsG3,k)*sp(k2,kp)))/3.;
                multTable[3][2]=(160*(sp(cepsG2,kp)*sp(epsG3,k) + sp(cepsG2,k)*sp(epsG3,kp) - sp(cepsG2,epsG3)*sp(k,kp)))/3.;
                multTable[3][3]=(160*(sp(cepsG3,kp)*sp(epsG3,k) + sp(cepsG3,k)*sp(epsG3,kp) - sp(cepsG3,epsG3)*sp(k,kp)))/3.;
                multTable[3][4]=(160*(sp(epsG3,kp)*sp(k,k3) - sp(epsG3,k3)*sp(k,kp) + sp(epsG3,k)*sp(k3,kp)))/3.;
                multTable[3][5]=0;
                multTable[3][6]=0;
                multTable[3][7]=0;
                multTable[3][8]=0;
                multTable[4][0]=(160*(sp(cepsPsi,kp)*sp(k,k3) - sp(cepsPsi,k3)*sp(k,kp) + sp(cepsPsi,k)*sp(k3,kp)))/3.;
                multTable[4][1]=(160*(-(sp(k,kp)*sp(k2,k3)) + sp(k,k3)*sp(k2,kp) + sp(k,k2)*sp(k3,kp)))/3.;
                multTable[4][2]=(160*(sp(cepsG2,kp)*sp(k,k3) - sp(cepsG2,k3)*sp(k,kp) + sp(cepsG2,k)*sp(k3,kp)))/3.;
                multTable[4][3]=(160*(sp(cepsG3,kp)*sp(k,k3) - sp(cepsG3,k3)*sp(k,kp) + sp(cepsG3,k)*sp(k3,kp)))/3.;
                multTable[4][4]=(320*sp(k,k3)*sp(k3,kp))/3.;
                multTable[4][5]=0;
                multTable[4][6]=0;
                multTable[4][7]=0;
                multTable[4][8]=0;
                multTable[5][0]=0;
                multTable[5][1]=0;
                multTable[5][2]=0;
                multTable[5][3]=0;
                multTable[5][4]=0;
                multTable[5][5]=192*sp(k,k2)*sp(k2,kp);
                multTable[5][6]=96*(-(sp(k,kp)*sp(k2,k3)) + sp(k,k3)*sp(k2,kp) + sp(k,k2)*sp(k3,kp));
                multTable[5][7]=96*(sp(cepsG3,kp)*sp(k,k2) - sp(cepsG3,k2)*sp(k,kp) + sp(cepsG3,k)*sp(k2,kp));
                multTable[5][8]=96*(sp(cepsG2,kp)*sp(k,k2) - sp(cepsG2,k2)*sp(k,kp) + sp(cepsG2,k)*sp(k2,kp));
                multTable[6][0]=0;
                multTable[6][1]=0;
                multTable[6][2]=0;
                multTable[6][3]=0;
                multTable[6][4]=0;
                multTable[6][5]=96*(-(sp(k,kp)*sp(k2,k3)) + sp(k,k3)*sp(k2,kp) + sp(k,k2)*sp(k3,kp));
                multTable[6][6]=192*sp(k,k3)*sp(k3,kp);
                multTable[6][7]=96*(sp(cepsG3,kp)*sp(k,k3) - sp(cepsG3,k3)*sp(k,kp) + sp(cepsG3,k)*sp(k3,kp));
                multTable[6][8]=96*(sp(cepsG2,kp)*sp(k,k3) - sp(cepsG2,k3)*sp(k,kp) + sp(cepsG2,k)*sp(k3,kp));
                multTable[7][0]=0;
                multTable[7][1]=0;
                multTable[7][2]=0;
                multTable[7][3]=0;
                multTable[7][4]=0;
                multTable[7][5]=96*(sp(epsG3,kp)*sp(k,k2) - sp(epsG3,k2)*sp(k,kp) + sp(epsG3,k)*sp(k2,kp));
                multTable[7][6]=96*(sp(epsG3,kp)*sp(k,k3) - sp(epsG3,k3)*sp(k,kp) + sp(epsG3,k)*sp(k3,kp));
                multTable[7][7]=96*(sp(cepsG3,kp)*sp(epsG3,k) + sp(cepsG3,k)*sp(epsG3,kp) - sp(cepsG3,epsG3)*sp(k,kp));
                multTable[7][8]=96*(sp(cepsG2,kp)*sp(epsG3,k) + sp(cepsG2,k)*sp(epsG3,kp) - sp(cepsG2,epsG3)*sp(k,kp));
                multTable[8][0]=0;
                multTable[8][1]=0;
                multTable[8][2]=0;
                multTable[8][3]=0;
                multTable[8][4]=0;
                multTable[8][5]=96*(sp(epsG2,kp)*sp(k,k2) - sp(epsG2,k2)*sp(k,kp) + sp(epsG2,k)*sp(k2,kp));
                multTable[8][6]=96*(sp(epsG2,kp)*sp(k,k3) - sp(epsG2,k3)*sp(k,kp) + sp(epsG2,k)*sp(k3,kp));
                multTable[8][7]=96*(sp(cepsG3,kp)*sp(epsG2,k) + sp(cepsG3,k)*sp(epsG2,kp) - sp(cepsG3,epsG2)*sp(k,kp));
                multTable[8][8]=96*(sp(cepsG2,kp)*sp(epsG2,k) + sp(cepsG2,k)*sp(epsG2,kp) - sp(cepsG2,epsG2)*sp(k,kp));

                 
                 for(int iAmp=0; iAmp<nAmps; ++iAmp) {
                     for(int ciAmp=0; ciAmp<nAmps; ++ciAmp) {
                         imatr2 += Amp[iAmp]*conj(Amp[ciAmp])*multTable[iAmp][ciAmp];
                     };
                 };                    
             }
         }
     }
//        for(int _i=0; _i<4; ++_i) epsG2[_i]=k2[_i];
     assert(is_zero(imatr2.imag()));
     assert(imatr2.real()>0);
     if(iEv<2) {
         cout<<"(*----- Debug print at iEv="<<iEv<<"-------*)"<<endl;
         println_4v_math("k",k);
         println_4v_math("kp",kp);
         println_4v_math("k1",k1);
         println_4v_math("k2",k2);
         println_4v_math("k3",k3);
         println_4v_math("pPsi",pPsi);
         println_4v_math("epsG2",epsG2);
         println_4v_math("epsG3",epsG3);
         println_4v_math("epsPsi",epsPsi);
         cout<<"MultTable"<<endl;
         for(int iA=0; iA<nAmps; ++iA) {
             for(int ciA=0; ciA<nAmps; ++ciA) {
                cout<<"\t {"<<iA<<","<<ciA<<","<<multTable[iA][ciA]<<"},"<<endl;
             }
         };
         cout<<"End of multtable"<<endl;
         cout<<"Amp"<<endl;
         for(int iA=0; iA<nAmps; ++iA) {
             cout<<"\t {"<<iA<<","<<Amp[iA]<<"},"<<endl;
         };
         cout<<"end of amp"<<endl;
         cout<<" imatr2="<<imatr2<<endl;
     }
     return imatr2.real();    
};



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

    int nBins=30;
    TFile file("out_pol_3S1_co.root","RECREATE");
    TNtuple tup("tup","tup","Q2:Y:pTpsi:W2:dsigma");
    TH1F hPT2("hPT2","hPT2",nBins,0,30); hPT2.Sumw2();
    TH1F hQ2("hQ2","hQ2",nBins,minQ2,maxQ2); hQ2.Sumw2();
    TH1F hZ("hZ","hZ",nBins,0,1); hZ.Sumw2();
    TH1F hY("hY","hY",nBins,0,1); hY.Sumw2();
    TH1F hW("hW","hW",nBins,sqrt(minW2),sqrt(maxW2)); hW.Sumw2();
    
    
    dbl_type Q2_scale;
    dbl_type P[4], k[4], kp[4], k1[4], k2[4], k3[4], pPsi[4];
    int nEv=1e3, nPassed=0, nNegative=0;
    dbl_type sum=0, dsigma, sigma=0;
    for(int iEv=0; iEv<nEv; ++iEv) {
        if( iEv % (nEv/100) == 0 && iEv>0) cout<<"---- Event "<<iEv<<" ("
                <<(int)(100.*iEv/nEv)<<" %) --- sigma="<<sigma*nEv/iEv<<" pb"<<endl;
        x=random->rand(minX,maxX);

        if(!kinematics(ramEP,k,kp,k1,k2,k3,pPsi)) continue;
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

        dbl_type matr2=getMatr2_pol(k,kp,k1,k2,k3,pPsi,nPassed);
        dbl_type wt=ramEP->wt;
        if(!(matr2>0)) {
            nNegative++;
            continue;
        };      

        dsigma=1;
        dsigma *= matr2*pdf*wt;
        dsigma *= 1./4/8; // polarizations and initial gluon spin
        dsigma *= 1./(2*x*pow(ecm,2)); // 4 I
        dsigma *= 1./nEv;
        dsigma *= picob;
        sigma += dsigma;

        tup.Fill(ramEP->Q2,ramEP->Y,pT(pPsi), W2, dsigma);
        hPT2.Fill(pT_squared(pPsi),dsigma);
        hQ2.Fill(ramEP->Q2,dsigma);
        hZ.Fill(z,dsigma);
        hY.Fill(ramEP->Y,dsigma);
        hW.Fill(sqrt(W2),dsigma);

        ++nPassed;
    };
    
    tup.Write(); hPT2.Write();
    file.Save();
    write_histogram_to_file(hPT2,"hPT2_pol_3S1_co.hst");
    write_histogram_to_file(hQ2,"hQ2_pol_3S1_co.hst");
    write_histogram_to_file(hZ,"hZ_pol_3S1_co.hst");
    write_histogram_to_file(hY,"hY_pol_3S1_co.hst");
    write_histogram_to_file(hW,"hW_pol_3S1_co.hst");

    cout<<" sigma="<<sigma<<" pb"<<endl;
    cout<<nPassed<<" ("<<(int)(100.*nPassed/nEv)<<"%) events passed"<<endl;
    cout<<nNegative<<" ("<<(int)(100.*nNegative/nEv)<<"%) events with negative matr2"<<endl;
    cout<<ramEP->nFault<<" faults in ramEP"<<endl;

 }