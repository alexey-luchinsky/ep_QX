#include "kinematics/RamboEP.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>

extern dbl_type mc, Mcc, ecm, PI, alpha, alphas, Opsi, x;
extern dbl_type P[4], k[4], kp[4], k1[4], k2[4], k3[4], pPsi[4];
cmplx epsG1[4], epsG2[4], epsG3[4], epsPsi[4], cepsG1[4], cepsG2[4], cepsG3[4], cepsPsi[4];
#define spp(k,kp,v1,v2) 8*4*(sp(k,v1)*sp(kp,v2)-sp(v1,v2)*sp(k,kp)+sp(kp,v1)*sp(k,v2))

const int max_nAmp=30;
cmplx multTable[max_nAmp][max_nAmp];
cmplx Amp[max_nAmp];



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


bool kinematics(RamboEP *ramEP) {
    if(!ramEP->next(kp,pPsi,k3,x)) return false;
    set_v4(k,ramEP->kIn);                // k
        assert(is_zero(mass2(k)));       
    set_v4(kp, ramEP->kOut);             // kp
        assert(is_zero(mass2(kp)));
    set_v4(k2,ramEP->Pg);                // k2
        assert(are_equal(sum_mass2(pPsi,k3),ramEP->W2));
    subtract(k,kp,k1);                  // k1=k-kp
        assert(are_equal(mass2(k1),-ramEP->Q2)); 
        assert(are_equal(mass2(pPsi),Mcc*Mcc)); 
        assert(is_zero(mass2(k3)));
        assert(are_equal(sum_mass2(k,k2),x*ecm*ecm));
        assert(are_equal(sum_mass2(kp,pPsi,k3),x*ecm*ecm));
    return true;
};


int calcMult_3S1_cs() {
    multTable[0][0]=spp(k,kp,epsPsi,cepsPsi);
    multTable[0][1]=spp(k,kp,epsPsi,k2);
    multTable[0][2]=spp(k,kp,epsPsi,cepsG2);
    multTable[0][3]=spp(k,kp,epsPsi,cepsG3);
    multTable[0][4]=spp(k,kp,epsPsi,k3);
    multTable[1][0]=spp(k,kp,k2,cepsPsi);
    multTable[1][1]=spp(k,kp,k2,k2);
    multTable[1][2]=spp(k,kp,k2,cepsG2);
    multTable[1][3]=spp(k,kp,k2,cepsG3);
    multTable[1][4]=spp(k,kp,k2,k3);
    multTable[2][0]=spp(k,kp,epsG2,cepsPsi);
    multTable[2][1]=spp(k,kp,epsG2,k2);
    multTable[2][2]=spp(k,kp,epsG2,cepsG2);
    multTable[2][3]=spp(k,kp,epsG2,cepsG3);
    multTable[2][4]=spp(k,kp,epsG2,k3);
    multTable[3][0]=spp(k,kp,epsG3,cepsPsi);
    multTable[3][1]=spp(k,kp,epsG3,k2);
    multTable[3][2]=spp(k,kp,epsG3,cepsG2);
    multTable[3][3]=spp(k,kp,epsG3,cepsG3);
    multTable[3][4]=spp(k,kp,epsG3,k3);
    multTable[4][0]=spp(k,kp,k3,cepsPsi);
    multTable[4][1]=spp(k,kp,k3,k2);
    multTable[4][2]=spp(k,kp,k3,cepsG2);
    multTable[4][3]=spp(k,kp,k3,cepsG3);
    multTable[4][4]=spp(k,kp,k3,k3);
    return 5;
}


void debug_print(int nAmps) {
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
}


int calcMult_3S1_co() {
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
    return 9;
}




int calcAmp_3S1_co() {
    dbl_type Q2=-mass2(k1);
    cmplx II=COMPLEX_ONE;
    Amp[0]=(-64*alpha*alphas*pow(PI,2)*(pow(sp(k2,k3),2)*sp(epsG2,epsG3) - sp(epsG2,pPsi)*sp(epsG3,k3)*sp(k2,k3) - sp(epsG2,pPsi)*sp(epsG3,pPsi)*sp(k2,k3) + sp(epsG2,epsG3)*sp(k2,k3)*sp(k2,pPsi) + sp(epsG2,k3)*((sp(epsG3,k3) + sp(epsG3,pPsi))*sp(k2,pPsi) - sp(epsG3,k2)*(sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))) + sp(epsG2,pPsi)*sp(epsG3,k2)*sp(k3,pPsi) - sp(epsG2,epsG3)*sp(k2,k3)*sp(k3,pPsi) - sp(epsG2,epsG3)*sp(k2,pPsi)*sp(k3,pPsi) + sp(epsG2,k2)*(sp(epsG3,k3)*sp(k2,k3) + sp(epsG3,pPsi)*sp(k2,k3) - sp(epsG3,k2)*sp(k3,pPsi)))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
    Amp[1]=(-16*alpha*alphas*pow(PI,2)*(4*sp(epsG2,k3)*sp(epsG3,epsPsi)*(sp(k2,k3) - sp(k3,pPsi)) + 4*sp(epsG2,epsPsi)*(sp(epsG3,k3)*sp(k2,k3) + sp(epsG3,pPsi)*sp(k2,k3) - sp(epsG3,k2)*sp(k3,pPsi)) + sp(epsG2,epsG3)*(-4*sp(epsPsi,k3)*(sp(k2,k3) - sp(k3,pPsi)) + sp(epsPsi,pPsi)*(-2*sp(k2,k3) + sp(k3,pPsi))))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
    Amp[2]=(-16*alpha*alphas*pow(PI,2)*(-4*pow(sp(k2,k3),2)*sp(epsG3,epsPsi) + 4*sp(epsG3,k2)*sp(epsPsi,k3)*sp(k2,k3) + 2*sp(epsG3,k2)*sp(epsPsi,pPsi)*sp(k2,k3) + sp(epsG3,k3)*(-4*sp(epsPsi,k2)*sp(k2,k3) + sp(epsPsi,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi))) + sp(epsG3,pPsi)*(-4*sp(epsPsi,k2)*sp(k2,k3) + sp(epsPsi,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi))) + 4*sp(epsG3,k2)*sp(epsPsi,k2)*sp(k3,pPsi) - 4*sp(epsG3,k2)*sp(epsPsi,k3)*sp(k3,pPsi) - 3*sp(epsG3,k2)*sp(epsPsi,pPsi)*sp(k3,pPsi) + 4*sp(epsG3,epsPsi)*sp(k2,k3)*sp(k3,pPsi))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
    Amp[3]=(16*alpha*alphas*pow(PI,2)*(4*pow(sp(k2,k3),2)*sp(epsG2,epsPsi) - 4*sp(epsG2,pPsi)*sp(epsPsi,k3)*sp(k2,k3) - 2*sp(epsG2,pPsi)*sp(epsPsi,pPsi)*sp(k2,k3) + 4*sp(epsG2,epsPsi)*sp(k2,k3)*sp(k2,pPsi) + sp(epsG2,k3)*(4*sp(epsPsi,k3)*sp(k2,pPsi) - 4*sp(epsPsi,k2)*(sp(k2,k3) + sp(k2,pPsi)) + sp(epsPsi,pPsi)*(2*sp(k2,k3) + 3*sp(k2,pPsi))) + sp(epsG2,k2)*(4*sp(epsPsi,k3)*sp(k2,k3) + sp(epsPsi,pPsi)*(2*sp(k2,k3) - sp(k3,pPsi))) + sp(epsG2,pPsi)*sp(epsPsi,pPsi)*sp(k3,pPsi))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
    Amp[4]=(-16*alpha*alphas*pow(PI,2)*(4*sp(epsG2,k2)*sp(epsG3,epsPsi)*sp(k2,k3) - 4*sp(epsG2,pPsi)*sp(epsG3,epsPsi)*sp(k2,k3) + 4*sp(epsG2,epsPsi)*sp(epsG3,k2)*sp(k2,k3) - 4*sp(epsG2,epsG3)*sp(epsPsi,k2)*sp(k2,k3) + 2*sp(epsG2,epsG3)*sp(epsPsi,pPsi)*sp(k2,k3) + 4*sp(epsG2,k3)*sp(epsG3,epsPsi)*sp(k2,pPsi) + 4*sp(epsG2,epsPsi)*sp(epsG3,k2)*sp(k2,pPsi) - 4*sp(epsG2,epsG3)*sp(epsPsi,k2)*sp(k2,pPsi) + sp(epsG2,epsG3)*sp(epsPsi,pPsi)*sp(k2,pPsi))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
    Amp[5]=(-16*alpha*alphas*II*pow(PI,2)*sp(epsG2,epsG3)*sp(epsPsi,pPsi)*(2*pow(sp(k2,k3),2) - sp(k2,k3)*sp(k3,pPsi) + sp(k2,pPsi)*sp(k3,pPsi))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,k3)*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
    Amp[6]=(-16*alpha*alphas*II*pow(PI,2)*sp(epsG2,epsG3)*sp(epsPsi,pPsi)*(2*pow(sp(k2,k3),2) + sp(k2,k3)*sp(k2,pPsi) + sp(k2,pPsi)*sp(k3,pPsi))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,k3)*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
    Amp[7]=(-16*alpha*alphas*II*pow(PI,2)*sp(epsPsi,pPsi)*(sp(epsG2,pPsi)*sp(k2,k3)*(-2*sp(k2,k3) + sp(k3,pPsi)) + sp(epsG2,k3)*(-2*pow(sp(k2,k3),2) + sp(k2,k3)*sp(k2,pPsi) - 2*sp(k2,pPsi)*sp(k3,pPsi)) + sp(epsG2,k2)*(2*pow(sp(k2,k3),2) - sp(k2,k3)*sp(k3,pPsi) + sp(k2,pPsi)*sp(k3,pPsi)))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,k3)*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
    Amp[8]=(-16*alpha*alphas*II*pow(PI,2)*sp(epsPsi,pPsi)*(sp(epsG3,pPsi)*sp(k2,k3)*(2*sp(k2,k3) + sp(k2,pPsi)) + sp(epsG3,k3)*(2*pow(sp(k2,k3),2) + sp(k2,k3)*sp(k2,pPsi) + sp(k2,pPsi)*sp(k3,pPsi)) - sp(epsG3,k2)*(2*pow(sp(k2,k3),2) + sp(k2,k3)*sp(k3,pPsi) + 2*sp(k2,pPsi)*sp(k3,pPsi)))*sqrt(mc*Opsi))/(3.*Q2*sp(k2,k3)*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi)*sqrt(3));
    return 9;
}

dbl_type getMatr2(std::function<int (void)> mtabl_func, std::function<int (void)> amp_func) {
    cmplx imatr2=COMPLEX_ZERO;
    int nAmps;
     for(int iPolPsi=0; iPolPsi<3; ++iPolPsi) {
         for(int iPolG2=0; iPolG2<2; ++iPolG2) {
             for(int iPolG3=0; iPolG3<2; ++iPolG3) {
                 set_gluon_polarization(iPolG2, k2, epsG2, cepsG2);
                 set_gluon_polarization(iPolG3, k3, epsG3, cepsG3);
                 set_psi_polarization(iPolPsi, pPsi, epsPsi, cepsPsi);
                 
                 nAmps=mtabl_func();
                 amp_func();
                 
                 for(int iAmp=0; iAmp<nAmps; ++iAmp) {
                     for(int ciAmp=0; ciAmp<nAmps; ++ciAmp) {
                         imatr2 += Amp[iAmp]*conj(Amp[ciAmp])*multTable[iAmp][ciAmp];
                     };
                 };                    
             }
         }
     }
     assert(is_zero(imatr2.imag()));
     assert(imatr2.real()>0);
     return imatr2.real();    
}

int calcAmp_3S1_cs() {
    dbl_type Q2=-mass2(k1);
    Amp[0]=(-128*alpha*alphas*mc*pow(PI,2)*(pow(sp(k2,k3),2)*sp(epsG2,epsG3) - sp(epsG2,pPsi)*sp(epsG3,k3)*sp(k2,k3) - sp(epsG2,pPsi)*sp(epsG3,pPsi)*sp(k2,k3) + sp(epsG2,epsG3)*sp(k2,k3)*sp(k2,pPsi) + sp(epsG2,k3)*((sp(epsG3,k3) + sp(epsG3,pPsi))*sp(k2,pPsi) - sp(epsG3,k2)*(sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))) + sp(epsG2,pPsi)*sp(epsG3,k2)*sp(k3,pPsi) - sp(epsG2,epsG3)*sp(k2,k3)*sp(k3,pPsi) - sp(epsG2,epsG3)*sp(k2,pPsi)*sp(k3,pPsi) + sp(epsG2,k2)*(sp(epsG3,k3)*sp(k2,k3) + sp(epsG3,pPsi)*sp(k2,k3) - sp(epsG3,k2)*sp(k3,pPsi)))*sqrt(Opsi/mc))/(9.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi));
    Amp[1]=(-32*alpha*alphas*mc*pow(PI,2)*(4*sp(epsG2,k3)*sp(epsG3,epsPsi)*(sp(k2,k3) - sp(k3,pPsi)) + 4*sp(epsG2,epsPsi)*(sp(epsG3,k3)*sp(k2,k3) + sp(epsG3,pPsi)*sp(k2,k3) - sp(epsG3,k2)*sp(k3,pPsi)) + sp(epsG2,epsG3)*(-4*sp(epsPsi,k3)*(sp(k2,k3) - sp(k3,pPsi)) + sp(epsPsi,pPsi)*(-2*sp(k2,k3) + sp(k3,pPsi))))*sqrt(Opsi/mc))/(9.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi));
    Amp[2]=(-32*alpha*alphas*mc*pow(PI,2)*(-4*pow(sp(k2,k3),2)*sp(epsG3,epsPsi) + 4*sp(epsG3,k2)*sp(epsPsi,k3)*sp(k2,k3) + 2*sp(epsG3,k2)*sp(epsPsi,pPsi)*sp(k2,k3) + sp(epsG3,k3)*(-4*sp(epsPsi,k2)*sp(k2,k3) + sp(epsPsi,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi))) + sp(epsG3,pPsi)*(-4*sp(epsPsi,k2)*sp(k2,k3) + sp(epsPsi,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi))) + 4*sp(epsG3,k2)*sp(epsPsi,k2)*sp(k3,pPsi) - 4*sp(epsG3,k2)*sp(epsPsi,k3)*sp(k3,pPsi) - 3*sp(epsG3,k2)*sp(epsPsi,pPsi)*sp(k3,pPsi) + 4*sp(epsG3,epsPsi)*sp(k2,k3)*sp(k3,pPsi))*sqrt(Opsi/mc))/(9.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi));
    Amp[3]=(32*alpha*alphas*mc*pow(PI,2)*(4*pow(sp(k2,k3),2)*sp(epsG2,epsPsi) - 4*sp(epsG2,pPsi)*sp(epsPsi,k3)*sp(k2,k3) - 2*sp(epsG2,pPsi)*sp(epsPsi,pPsi)*sp(k2,k3) + 4*sp(epsG2,epsPsi)*sp(k2,k3)*sp(k2,pPsi) + sp(epsG2,k3)*(4*sp(epsPsi,k3)*sp(k2,pPsi) - 4*sp(epsPsi,k2)*(sp(k2,k3) + sp(k2,pPsi)) + sp(epsPsi,pPsi)*(2*sp(k2,k3) + 3*sp(k2,pPsi))) + sp(epsG2,k2)*(4*sp(epsPsi,k3)*sp(k2,k3) + sp(epsPsi,pPsi)*(2*sp(k2,k3) - sp(k3,pPsi))) + sp(epsG2,pPsi)*sp(epsPsi,pPsi)*sp(k3,pPsi))*sqrt(Opsi/mc))/(9.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi));
    Amp[4]=(-32*alpha*alphas*mc*pow(PI,2)*(4*sp(epsG2,k2)*sp(epsG3,epsPsi)*sp(k2,k3) - 4*sp(epsG2,pPsi)*sp(epsG3,epsPsi)*sp(k2,k3) + 4*sp(epsG2,epsPsi)*sp(epsG3,k2)*sp(k2,k3) - 4*sp(epsG2,epsG3)*sp(epsPsi,k2)*sp(k2,k3) + 2*sp(epsG2,epsG3)*sp(epsPsi,pPsi)*sp(k2,k3) + 4*sp(epsG2,k3)*sp(epsG3,epsPsi)*sp(k2,pPsi) + 4*sp(epsG2,epsPsi)*sp(epsG3,k2)*sp(k2,pPsi) - 4*sp(epsG2,epsG3)*sp(epsPsi,k2)*sp(k2,pPsi) + sp(epsG2,epsG3)*sp(epsPsi,pPsi)*sp(k2,pPsi))*sqrt(Opsi/mc))/(9.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi));
    return 5;
};

/*
dbl_type getMatr2_3S1_cs(dbl_type(&k)[4], dbl_type(&kp)[4], dbl_type(&k1)[4], dbl_type(&k2)[4], 
                dbl_type(&k3)[4], dbl_type(&pPsi)[4], int iEv) {
    return getMatr2(calcMult_3S1_cs, calcAmp_3S1_cs);
};
*/

dbl_type getMatr2_3S1_cs(dbl_type(&k)[4], dbl_type(&kp)[4], dbl_type(&k1)[4], dbl_type(&k2)[4], 
                dbl_type(&k3)[4], dbl_type(&pPsi)[4], int iEv) {
    dbl_type matr2=0;
    dbl_type Q2=-mass2(k1);
    for(int iG1=0; iG1<2; ++iG1) {
        lepton_current(iG1,k,kp,epsG1);
//        println_4v("epsG1",epsG1);
        for(int iG2=0; iG2<2; ++iG2) {
            set_gluon_polarization(iG2, k2, epsG2);
            for(int iG3=0; iG3<2; ++iG3) {
                set_gluon_polarization(iG3,k3,epsG3);
                for(int iPsi=0; iPsi<3; ++iPsi) {
                    set_psi_polarization(iPsi,pPsi, epsPsi);
                    cmplx amp=(-64*alpha*alphas*mc*pow(PI,2)*(4*pow(sp(k2,k3),2)*sp(epsG1,epsPsi)*sp(epsG2,epsG3) - 4*pow(sp(k2,k3),2)*sp(epsG1,epsG3)*sp(epsG2,epsPsi) - 4*pow(sp(k2,k3),2)*sp(epsG1,epsG2)*sp(epsG3,epsPsi) - 4*sp(epsG1,epsPsi)*sp(epsG2,k3)*sp(epsG3,k2)*sp(k2,k3) + 4*sp(epsG1,epsPsi)*sp(epsG2,k2)*sp(epsG3,k3)*sp(k2,k3) - 4*sp(epsG1,epsPsi)*sp(epsG2,pPsi)*sp(epsG3,k3)*sp(k2,k3) + 4*sp(epsG1,epsPsi)*sp(epsG2,k2)*sp(epsG3,pPsi)*sp(k2,k3) - 4*sp(epsG1,epsPsi)*sp(epsG2,pPsi)*sp(epsG3,pPsi)*sp(k2,k3) + 4*sp(epsG1,epsG3)*sp(epsG2,k3)*sp(epsPsi,k2)*sp(k2,k3) - 4*sp(epsG1,epsG2)*sp(epsG3,k3)*sp(epsPsi,k2)*sp(k2,k3) - 4*sp(epsG1,epsG2)*sp(epsG3,pPsi)*sp(epsPsi,k2)*sp(k2,k3) - 4*sp(epsG1,epsG3)*sp(epsG2,k2)*sp(epsPsi,k3)*sp(k2,k3) + 4*sp(epsG1,epsG3)*sp(epsG2,pPsi)*sp(epsPsi,k3)*sp(k2,k3) + 4*sp(epsG1,epsG2)*sp(epsG3,k2)*sp(epsPsi,k3)*sp(k2,k3) - 2*sp(epsG1,epsG3)*sp(epsG2,k2)*sp(epsPsi,pPsi)*sp(k2,k3) - 2*sp(epsG1,epsG3)*sp(epsG2,k3)*sp(epsPsi,pPsi)*sp(k2,k3) + 2*sp(epsG1,epsG3)*sp(epsG2,pPsi)*sp(epsPsi,pPsi)*sp(k2,k3) + 2*sp(epsG1,epsG2)*sp(epsG3,k2)*sp(epsPsi,pPsi)*sp(k2,k3) + 2*sp(epsG1,epsG2)*sp(epsG3,k3)*sp(epsPsi,pPsi)*sp(k2,k3) + 2*sp(epsG1,epsG2)*sp(epsG3,pPsi)*sp(epsPsi,pPsi)*sp(k2,k3) - 4*sp(epsG1,epsPsi)*sp(epsG2,k3)*sp(epsG3,k2)*sp(k2,pPsi) + 4*sp(epsG1,epsPsi)*sp(epsG2,k3)*sp(epsG3,k3)*sp(k2,pPsi) + 4*sp(epsG1,epsPsi)*sp(epsG2,k3)*sp(epsG3,pPsi)*sp(k2,pPsi) + 4*sp(epsG1,epsG3)*sp(epsG2,k3)*sp(epsPsi,k2)*sp(k2,pPsi) - 4*sp(epsG1,epsG3)*sp(epsG2,k3)*sp(epsPsi,k3)*sp(k2,pPsi) - 3*sp(epsG1,epsG3)*sp(epsG2,k3)*sp(epsPsi,pPsi)*sp(k2,pPsi) + sp(epsG1,epsG2)*sp(epsG3,k3)*sp(epsPsi,pPsi)*sp(k2,pPsi) + sp(epsG1,epsG2)*sp(epsG3,pPsi)*sp(epsPsi,pPsi)*sp(k2,pPsi) + 4*sp(epsG1,epsPsi)*sp(epsG2,epsG3)*sp(k2,k3)*sp(k2,pPsi) - 4*sp(epsG1,epsG3)*sp(epsG2,epsPsi)*sp(k2,k3)*sp(k2,pPsi) + sp(epsG1,k3)*(4*sp(epsG2,k2)*sp(epsG3,epsPsi)*sp(k2,k3) - 4*sp(epsG2,pPsi)*sp(epsG3,epsPsi)*sp(k2,k3) + 4*sp(epsG2,epsPsi)*sp(epsG3,k2)*sp(k2,k3) - 4*sp(epsG2,epsG3)*sp(epsPsi,k2)*sp(k2,k3) + 2*sp(epsG2,epsG3)*sp(epsPsi,pPsi)*sp(k2,k3) + 4*sp(epsG2,k3)*sp(epsG3,epsPsi)*sp(k2,pPsi) + 4*sp(epsG2,epsPsi)*sp(epsG3,k2)*sp(k2,pPsi) - 4*sp(epsG2,epsG3)*sp(epsPsi,k2)*sp(k2,pPsi) + sp(epsG2,epsG3)*sp(epsPsi,pPsi)*sp(k2,pPsi)) - 4*sp(epsG1,epsPsi)*sp(epsG2,k2)*sp(epsG3,k2)*sp(k3,pPsi) + 4*sp(epsG1,epsPsi)*sp(epsG2,k3)*sp(epsG3,k2)*sp(k3,pPsi) + 4*sp(epsG1,epsPsi)*sp(epsG2,pPsi)*sp(epsG3,k2)*sp(k3,pPsi) + 4*sp(epsG1,epsG2)*sp(epsG3,k2)*sp(epsPsi,k2)*sp(k3,pPsi) - 4*sp(epsG1,epsG2)*sp(epsG3,k2)*sp(epsPsi,k3)*sp(k3,pPsi) + sp(epsG1,epsG3)*sp(epsG2,k2)*sp(epsPsi,pPsi)*sp(k3,pPsi) - sp(epsG1,epsG3)*sp(epsG2,pPsi)*sp(epsPsi,pPsi)*sp(k3,pPsi) - 3*sp(epsG1,epsG2)*sp(epsG3,k2)*sp(epsPsi,pPsi)*sp(k3,pPsi) - 4*sp(epsG1,epsPsi)*sp(epsG2,epsG3)*sp(k2,k3)*sp(k3,pPsi) + 4*sp(epsG1,epsG2)*sp(epsG3,epsPsi)*sp(k2,k3)*sp(k3,pPsi) - 4*sp(epsG1,epsPsi)*sp(epsG2,epsG3)*sp(k2,pPsi)*sp(k3,pPsi) + sp(epsG1,k2)*(4*sp(epsG2,k3)*sp(epsG3,epsPsi)*(sp(k2,k3) - sp(k3,pPsi)) + 4*sp(epsG2,epsPsi)*(sp(epsG3,k3)*sp(k2,k3) + sp(epsG3,pPsi)*sp(k2,k3) - sp(epsG3,k2)*sp(k3,pPsi)) + sp(epsG2,epsG3)*(-4*sp(epsPsi,k3)*(sp(k2,k3) - sp(k3,pPsi)) + sp(epsPsi,pPsi)*(-2*sp(k2,k3) + sp(k3,pPsi)))))*sqrt(2)*sqrt(Opsi/mc))/(9.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi));
                    matr2 += pow(abs(amp),2);
//                    cout<<"matr2="<<matr2<<endl;
                }
            }
        }
    }
    return matr2;
}

dbl_type getMatr2_3S1_co(dbl_type(&k)[4], dbl_type(&kp)[4], dbl_type(&k1)[4], dbl_type(&k2)[4], 
                dbl_type(&k3)[4], dbl_type(&pPsi)[4], int iEv) {
    return getMatr2(calcMult_3S1_co, calcAmp_3S1_co);
};

