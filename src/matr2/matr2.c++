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
extern bool gauge;
cmplx epsG1[4], epsG2[4], epsG3[4], epsPsi[4];




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




//void debug_print(int nAmps) {
//    println_4v_math("k",k);
//    println_4v_math("kp",kp);
//    println_4v_math("k1",k1);
//    println_4v_math("k2",k2);
//    println_4v_math("k3",k3);
//    println_4v_math("pPsi",pPsi);
//    println_4v_math("epsG2",epsG2);
//    println_4v_math("epsG3",epsG3);
//    println_4v_math("epsPsi",epsPsi);
//    cout<<"MultTable"<<endl;
//    for(int iA=0; iA<nAmps; ++iA) {
//        for(int ciA=0; ciA<nAmps; ++ciA) {
//           cout<<"\t {"<<iA<<","<<ciA<<","<<multTable[iA][ciA]<<"},"<<endl;
//        }
//    };
//    cout<<"End of multtable"<<endl;
//    cout<<"Amp"<<endl;
//    for(int iA=0; iA<nAmps; ++iA) {
//        cout<<"\t {"<<iA<<","<<Amp[iA]<<"},"<<endl;
//    };
//    cout<<"end of amp"<<endl;
//}

dbl_type getMatr2_3S1_cs(dbl_type(&k)[4], dbl_type(&kp)[4], dbl_type(&k1)[4], dbl_type(&k2)[4], 
                dbl_type(&k3)[4], dbl_type(&pPsi)[4], int iEv) {
    dbl_type matr2=0;
    dbl_type Q2=-mass2(k1);
    for(int iG1=0; iG1<2; ++iG1) {
        lepton_current(iG1,k,kp,epsG1);
//        println_4v("epsG1",epsG1);
        for(int iG2=0; iG2<2; ++iG2) {
            set_gluon_polarization(iG2, k2, epsG2);
            if(gauge) for(int i=0; i<4; ++i) epsG2[i]=k2[i];
            for(int iG3=0; iG3<2; ++iG3) {
                set_gluon_polarization(iG3,k3,epsG3);
                for(int iPsi=0; iPsi<3; ++iPsi) {
                    set_psi_polarization(iPsi,pPsi, epsPsi);
                    cmplx amp=COMPLEX_ZERO;
                    amp += sp(epsG1,epsPsi)*(4*(sp(epsG2,k3)*(sp(epsG3,pPsi)*sp(k2,pPsi) - sp(epsG3,k2)*(sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))) + sp(epsG2,epsG3)*(sp(k2,k3) + sp(k2,pPsi))*(sp(k2,k3) - sp(k3,pPsi)) + sp(epsG2,pPsi)*(-(sp(epsG3,pPsi)*sp(k2,k3)) + sp(epsG3,k2)*sp(k3,pPsi))));
                    amp += sp(epsG1,epsG3)*(4*(-(sp(k2,k3)*(-(sp(epsG2,pPsi)*sp(epsPsi,k3)) + sp(epsG2,epsPsi)*(sp(k2,k3) + sp(k2,pPsi)))) + sp(epsG2,k3)*(-(sp(epsPsi,k3)*sp(k2,pPsi)) + sp(epsPsi,k2)*(sp(k2,k3) + sp(k2,pPsi)))));
                    amp += sp(epsG1,epsG2)*(4*(-(sp(epsG3,pPsi)*sp(epsPsi,k2)*sp(k2,k3)) + sp(epsG3,epsPsi)*sp(k2,k3)*(-sp(k2,k3) + sp(k3,pPsi)) + sp(epsG3,k2)*(sp(epsPsi,k3)*(sp(k2,k3) - sp(k3,pPsi)) + sp(epsPsi,k2)*sp(k3,pPsi))));
                    amp += sp(epsG1,k3)*(-4*(sp(epsG2,pPsi)*sp(epsG3,epsPsi)*sp(k2,k3) + sp(epsG2,epsG3)*sp(epsPsi,k2)*sp(k2,k3) - sp(epsG2,k3)*sp(epsG3,epsPsi)*sp(k2,pPsi) + sp(epsG2,epsG3)*sp(epsPsi,k2)*sp(k2,pPsi) - sp(epsG2,epsPsi)*sp(epsG3,k2)*(sp(k2,k3) + sp(k2,pPsi))));
                    amp += sp(epsG1,k2)*(4*(sp(epsG2,k3)*sp(epsG3,epsPsi)*(sp(k2,k3) - sp(k3,pPsi)) + sp(epsG2,epsG3)*sp(epsPsi,k3)*(-sp(k2,k3) + sp(k3,pPsi)) + sp(epsG2,epsPsi)*(sp(epsG3,pPsi)*sp(k2,k3) - sp(epsG3,k2)*sp(k3,pPsi))));
                    amp *= (-64*alpha*alphas*mc*pow(PI,2)*sqrt(2)*sqrt(Opsi/mc))/(9.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi));
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
    dbl_type matr2=0;
    dbl_type Q2=-mass2(k1);
    for(int iG1=0; iG1<2; ++iG1) {
        lepton_current(iG1,k,kp,epsG1);
//        println_4v("epsG1",epsG1);
        for(int iG2=0; iG2<2; ++iG2) {
            set_gluon_polarization(iG2, k2, epsG2);
            if(gauge) for(int i=0; i<4; ++i) epsG2[i]=k2[i];
            for(int iG3=0; iG3<2; ++iG3) {
                set_gluon_polarization(iG3,k3,epsG3);
                for(int iPsi=0; iPsi<3; ++iPsi) {
                    set_psi_polarization(iPsi,pPsi, epsPsi);
                        cmplx sund=COMPLEX_ZERO;
                        sund += sp(epsG1,epsPsi)*(sp(epsG2,k3)*(sp(epsG3,pPsi)*sp(k2,pPsi) - sp(epsG3,k2)*(sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))) + sp(epsG2,epsG3)*(sp(k2,k3) + sp(k2,pPsi))*(sp(k2,k3) - sp(k3,pPsi)) + sp(epsG2,pPsi)*(-(sp(epsG3,pPsi)*sp(k2,k3)) + sp(epsG3,k2)*sp(k3,pPsi)));
                        sund += sp(epsG1,epsG3)*(-(sp(k2,k3)*(-(sp(epsG2,pPsi)*sp(epsPsi,k3)) + sp(epsG2,epsPsi)*(sp(k2,k3) + sp(k2,pPsi)))) + sp(epsG2,k3)*(-(sp(epsPsi,k3)*sp(k2,pPsi)) + sp(epsPsi,k2)*(sp(k2,k3) + sp(k2,pPsi))));
                        sund += sp(epsG1,epsG2)*(-(sp(epsG3,pPsi)*sp(epsPsi,k2)*sp(k2,k3)) + sp(epsG3,epsPsi)*sp(k2,k3)*(-sp(k2,k3) + sp(k3,pPsi)) + sp(epsG3,k2)*(sp(epsPsi,k3)*(sp(k2,k3) - sp(k3,pPsi)) + sp(epsPsi,k2)*sp(k3,pPsi)));
                        sund += sp(epsG1,k3)*(-(sp(epsG2,pPsi)*sp(epsG3,epsPsi)*sp(k2,k3)) - sp(epsG2,epsG3)*sp(epsPsi,k2)*sp(k2,k3) + sp(epsG2,k3)*sp(epsG3,epsPsi)*sp(k2,pPsi) - sp(epsG2,epsG3)*sp(epsPsi,k2)*sp(k2,pPsi) + sp(epsG2,epsPsi)*sp(epsG3,k2)*(sp(k2,k3) + sp(k2,pPsi)));
                        sund += sp(epsG1,k2)*(sp(epsG2,k3)*sp(epsG3,epsPsi)*(sp(k2,k3) - sp(k3,pPsi)) + sp(epsG2,epsG3)*sp(epsPsi,k3)*(-sp(k2,k3) + sp(k3,pPsi)) + sp(epsG2,epsPsi)*(sp(epsG3,pPsi)*sp(k2,k3) - sp(epsG3,k2)*sp(k3,pPsi)));
                        sund *= (-64*alpha*alphas*mc*pow(PI,2)*sqrt(Opsi/mc))/(9.*Q2*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi));
                        sund *= sqrt(1.*40/3);  // SUND^2 
                        matr2 += pow(abs(sund),2);
               }
            }
        }
    }
    return matr2;
}

dbl_type getMatr2_1S0_co(dbl_type(&k)[4], dbl_type(&kp)[4], dbl_type(&k1)[4], dbl_type(&k2)[4], 
                dbl_type(&k3)[4], dbl_type(&pPsi)[4], int iEv) {
    dbl_type matr2=0;
    dbl_type Q2=-mass2(k1);
    for(int iG1=0; iG1<2; ++iG1) {
        lepton_current(iG1,k,kp,epsG1);
//        println_4v("epsG1",epsG1);
        for(int iG2=0; iG2<2; ++iG2) {
            set_gluon_polarization(iG2, k2, epsG2);
            if(gauge) for(int i=0; i<4; ++i) epsG2[i]=k2[i];
            for(int iG3=0; iG3<2; ++iG3) {
                set_gluon_polarization(iG3,k3,epsG3);
                cmplx sunf=COMPLEX_ZERO;
                sunf += lcv(epsG1,k2,k3,pPsi)*(sp(epsG2,epsG3)*(2*pow(sp(k2,k3),2) - sp(k2,k3)*sp(k2,pPsi) + 2*sp(k2,pPsi)*sp(k3,pPsi)));
                sunf += lcv(epsG1,epsG3,k2,k3)*(2*sp(epsG2,pPsi)*sp(k2,k3)*(-sp(k2,k3) + sp(k3,pPsi)));
                sunf += lcv(epsG1,epsG3,k3,pPsi)*(-2*(sp(epsG2,pPsi)*sp(k2,k3) - sp(epsG2,k3)*sp(k2,pPsi))*(sp(k2,k3) - sp(k3,pPsi)));
                sunf += lcv(epsG1,epsG2,k2,k3)*(-(sp(epsG3,pPsi)*sp(k2,k3)*(2*sp(k2,k3) + sp(k2,pPsi))));
                sunf += lcv(epsG1,epsG2,k2,pPsi)*(-2*(sp(k2,k3) + sp(k2,pPsi))*(sp(epsG3,pPsi)*sp(k2,k3) - sp(epsG3,k2)*sp(k3,pPsi)));
                sunf += lcv(epsG1,epsG2,epsG3,pPsi)*(-2*sp(k2,k3)*(sp(k2,k3) + sp(k2,pPsi))*(sp(k2,k3) - sp(k3,pPsi)));
                sunf += sp(epsG1,epsG3)*(lcv(epsG2,k2,k3,pPsi)*sp(k2,k3)*sp(k2,pPsi));
                sunf += lcv(epsG1,epsG3,k2,pPsi)*(sp(epsG2,k3)*(2*sp(k2,pPsi)*sp(k3,pPsi) + sp(k2,k3)*(-sp(k2,pPsi) + sp(k3,pPsi))));
                sunf += lcv(epsG1,epsG2,k3,pPsi)*(sp(epsG3,k2)*(sp(k2,k3)*(2*sp(k2,pPsi) - sp(k3,pPsi)) - 2*sp(k2,pPsi)*sp(k3,pPsi)));
                sunf += lcv(epsG1,epsG2,epsG3,k3)*(sp(k2,k3)*sp(k2,pPsi)*(sp(k2,pPsi) + sp(k3,pPsi)));
                sunf += sp(epsG1,k2)*(lcv(epsG2,epsG3,k3,pPsi)*sp(k2,k3)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi)));
                sunf += sp(epsG1,k3)*(lcv(epsG2,epsG3,k2,pPsi)*sp(k2,k3)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi)));
                sunf += lcv(epsG1,epsG2,epsG3,k2)*(-(sp(k2,k3)*sp(k3,pPsi)*(sp(k2,pPsi) + sp(k3,pPsi))));
                sunf *= (-16*alpha*alphas*pow(PI,2)*sqrt(Opsi/mc))/(9.*Q2*sp(k2,k3)*sp(k2,pPsi)*(2*sp(k2,k3) + sp(k2,pPsi) - sp(k3,pPsi))*sp(k3,pPsi));
                sunf *= sqrt(24.);
                matr2 += pow(abs(sunf),2);
//                    cout<<"matr2="<<matr2<<endl;
            }
        }
    }
    return matr2;
}


