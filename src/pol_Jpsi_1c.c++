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

const int nAmps=5;
cmplx Amp[nAmps];
cmplx multTable[nAmps][nAmps];
cmplx epsG2[4], epsG3[4], epsPsi[4];

void calcAmps(dbl_type(&k)[4], dbl_type(&kp)[4], dbl_type(&k1)[4], dbl_type(&k2)[4], 
                dbl_type(&k3)[4], dbl_type(&p)[4]) {
    dbl_type Q2=-mass2(k1);
	 Amp[0]=(-64*alpha*alphas*mc*pow(PI,2)*(-4*pow(sp(k2,k3),2)*sp(epsG3,epsPsi) + 4*sp(epsG3,k2)*sp(epsPsi,k3)*sp(k2,k3) + 2*sp(epsG3,k2)*sp(epsPsi,p)*sp(k2,k3) + sp(epsG3,k3)*(-4*sp(epsPsi,k2)*sp(k2,k3) + sp(epsPsi,p)*(2*sp(k2,k3) + sp(k2,p))) + sp(epsG3,p)*(-4*sp(epsPsi,k2)*sp(k2,k3) + sp(epsPsi,p)*(2*sp(k2,k3) + sp(k2,p))) + 4*sp(epsG3,k2)*sp(epsPsi,k2)*sp(k3,p) - 4*sp(epsG3,k2)*sp(epsPsi,k3)*sp(k3,p) - 3*sp(epsG3,k2)*sp(epsPsi,p)*sp(k3,p) + 4*sp(epsG3,epsPsi)*sp(k2,k3)*sp(k3,p))*sqrt(2)*sqrt(Opsi/mc))/(9.*Q2*sp(k2,p)*(2*sp(k2,k3) + sp(k2,p) - sp(k3,p))*sp(k3,p));
	 Amp[1]=(64*alpha*alphas*mc*pow(PI,2)*(4*pow(sp(k2,k3),2)*sp(epsG2,epsPsi) - 4*sp(epsG2,p)*sp(epsPsi,k3)*sp(k2,k3) - 2*sp(epsG2,p)*sp(epsPsi,p)*sp(k2,k3) + 4*sp(epsG2,epsPsi)*sp(k2,k3)*sp(k2,p) + sp(epsG2,k3)*(4*sp(epsPsi,k3)*sp(k2,p) - 4*sp(epsPsi,k2)*(sp(k2,k3) + sp(k2,p)) + sp(epsPsi,p)*(2*sp(k2,k3) + 3*sp(k2,p))) + sp(epsG2,k2)*(4*sp(epsPsi,k3)*sp(k2,k3) + sp(epsPsi,p)*(2*sp(k2,k3) - sp(k3,p))) + sp(epsG2,p)*sp(epsPsi,p)*sp(k3,p))*sqrt(2)*sqrt(Opsi/mc))/(9.*Q2*sp(k2,p)*(2*sp(k2,k3) + sp(k2,p) - sp(k3,p))*sp(k3,p));
	 Amp[2]=(-256*alpha*alphas*mc*pow(PI,2)*(pow(sp(k2,k3),2)*sp(epsG2,epsG3) - sp(epsG2,p)*sp(epsG3,k3)*sp(k2,k3) - sp(epsG2,p)*sp(epsG3,p)*sp(k2,k3) + sp(epsG2,epsG3)*sp(k2,k3)*sp(k2,p) + sp(epsG2,k3)*((sp(epsG3,k3) + sp(epsG3,p))*sp(k2,p) - sp(epsG3,k2)*(sp(k2,k3) + sp(k2,p) - sp(k3,p))) + sp(epsG2,p)*sp(epsG3,k2)*sp(k3,p) - sp(epsG2,epsG3)*sp(k2,k3)*sp(k3,p) - sp(epsG2,epsG3)*sp(k2,p)*sp(k3,p) + sp(epsG2,k2)*(sp(epsG3,k3)*sp(k2,k3) + sp(epsG3,p)*sp(k2,k3) - sp(epsG3,k2)*sp(k3,p)))*sqrt(2)*sqrt(Opsi/mc))/(9.*Q2*sp(k2,p)*(2*sp(k2,k3) + sp(k2,p) - sp(k3,p))*sp(k3,p));
	 Amp[3]=(-64*alpha*alphas*mc*pow(PI,2)*(4*sp(epsG2,k3)*sp(epsG3,epsPsi)*(sp(k2,k3) - sp(k3,p)) + 4*sp(epsG2,epsPsi)*(sp(epsG3,k3)*sp(k2,k3) + sp(epsG3,p)*sp(k2,k3) - sp(epsG3,k2)*sp(k3,p)) + sp(epsG2,epsG3)*(-4*sp(epsPsi,k3)*(sp(k2,k3) - sp(k3,p)) + sp(epsPsi,p)*(-2*sp(k2,k3) + sp(k3,p))))*sqrt(2)*sqrt(Opsi/mc))/(9.*Q2*sp(k2,p)*(2*sp(k2,k3) + sp(k2,p) - sp(k3,p))*sp(k3,p));
	 Amp[4]=(-64*alpha*alphas*mc*pow(PI,2)*(4*sp(epsG2,k2)*sp(epsG3,epsPsi)*sp(k2,k3) - 4*sp(epsG2,p)*sp(epsG3,epsPsi)*sp(k2,k3) + 4*sp(epsG2,epsPsi)*sp(epsG3,k2)*sp(k2,k3) - 4*sp(epsG2,epsG3)*sp(epsPsi,k2)*sp(k2,k3) + 2*sp(epsG2,epsG3)*sp(epsPsi,p)*sp(k2,k3) + 4*sp(epsG2,k3)*sp(epsG3,epsPsi)*sp(k2,p) + 4*sp(epsG2,epsPsi)*sp(epsG3,k2)*sp(k2,p) - 4*sp(epsG2,epsG3)*sp(epsPsi,k2)*sp(k2,p) + sp(epsG2,epsG3)*sp(epsPsi,p)*sp(k2,p))*sqrt(2)*sqrt(Opsi/mc))/(9.*Q2*sp(k2,p)*(2*sp(k2,k3) + sp(k2,p) - sp(k3,p))*sp(k3,p));
}

cmplx spp(dbl_type (&k)[4], dbl_type (&kp)[4], cmplx (&v1)[4], cmplx (&v2)[4]) {
    cmplx cv2[4];
    for(int i=0; i<3; ++i) cv2[i]=conj(v2[i]);
    return 4*(sp(k,v1)*sp(kp,cv2)-sp(v1,cv2)*sp(k,kp)+sp(kp,v1)*sp(k,cv2));
}

cmplx spp(dbl_type (&k)[4], dbl_type (&kp)[4], dbl_type (&v1)[4], cmplx (&v2)[4]) {
    cmplx cv2[4];
    for(int i=0; i<3; ++i) cv2[i]=conj(v2[i]);
    return 4*(sp(k,v1)*sp(kp,cv2)-sp(v1,cv2)*sp(k,kp)+sp(kp,v1)*sp(k,cv2));
}

cmplx spp(dbl_type (&k)[4], dbl_type (&kp)[4], cmplx (&v1)[4], dbl_type (&v2)[4]) {
    return 4*(sp(k,v1)*sp(kp,v2)-sp(v1,v2)*sp(k,kp)+sp(kp,v1)*sp(k,v2));
}

cmplx spp(dbl_type (&k)[4], dbl_type (&kp)[4], dbl_type (&v1)[4], dbl_type (&v2)[4]) {
    return 4*(sp(k,v1)*sp(kp,v2)-sp(v1,v2)*sp(k,kp)+sp(kp,v1)*sp(k,v2));
}

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


int main(void) {
    string pdf_set_name="cteq6l1"; cout<<" pdfse="<<pdf_set_name<<endl;
    dbl_type xi=1; cout<<" xi="<<xi<<endl;
    lhapdf_pdf=LHAPDF::mkPDF(pdf_set_name,0);
    Random *random=new Random();

    // experimental conditions from Kinehl
    ecm=300;
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
    
    

    dbl_type Q2_scale;
    dbl_type P[4], k[4], kp[4], k1[4], k2[4], k3[4], pPsi[4];
    int nEv=10, nPassed=0, nNegative=0;
    dbl_type sum=0, dsigma, sigma=0;
    for(int iEv=0; iEv<nEv; ++iEv) {
        if( iEv % (nEv/10) == 0 && iEv>0) cout<<"---- Event "<<iEv<<" ("
                <<(int)(100.*iEv/nEv)<<" %) --- sigma="<<sigma*nEv/iEv<<" pb"<<endl;
        x=random->rand(minX,maxX);

        if(!kinematics(ramEP,k,kp,k1,k2,k3,pPsi)) continue;
        ramEP->wt *= maxX-minX;
        dbl_type Q2=-mass2(k1);
        set_v4(P,ramEP->P);
        dbl_type Y=sp(k1,P)/sp(k,P);
        dbl_type z=sp(pPsi,P)/sp(k1,P);
        dbl_type W2=sum_mass2(k1,P);

        dbl_type matr2=0;
        for(int iPolPsi=0; iPolPsi<3; ++iPolPsi) {
            for(int iPolG2=0; iPolG2<2; ++iPolG2) {
                for(int iPolG3=0; iPolG3<2; ++iPolG3) {
                    set_gluon_polarization(iPolG2, k2, epsG2);
//                    for(int _i=0; _i<4; ++_i) epsG2[_i]=k2[_i];
                    set_gluon_polarization(iPolG3, k3, epsG3);
                    set_psi_polarization(iPolPsi, pPsi, epsPsi);
                    calcAmps(k,kp,k1,k2,k3,pPsi);
                    multTable[0][0]=spp(k,kp,epsG2,epsG2);
                    multTable[0][1]=spp(k,kp,epsG2,epsG3);
                    multTable[0][2]=spp(k,kp,epsG2,epsPsi);
                    multTable[0][3]=spp(k,kp,epsG2,k2);
                    multTable[0][4]=spp(k,kp,epsG2,k3);
                    multTable[1][0]=spp(k,kp,epsG3,epsG2);
                    multTable[1][1]=spp(k,kp,epsG3,epsG3);
                    multTable[1][2]=spp(k,kp,epsG3,epsPsi);
                    multTable[1][3]=spp(k,kp,epsG3,k2);
                    multTable[1][4]=spp(k,kp,epsG3,k3);
                    multTable[2][0]=spp(k,kp,epsPsi,epsG2);
                    multTable[2][1]=spp(k,kp,epsPsi,epsG3);
                    multTable[2][2]=spp(k,kp,epsPsi,epsPsi);
                    multTable[2][3]=spp(k,kp,epsPsi,k2);
                    multTable[2][4]=spp(k,kp,epsPsi,k3);
                    multTable[3][0]=spp(k,kp,k2,epsG2);
                    multTable[3][1]=spp(k,kp,k2,epsG3);
                    multTable[3][2]=spp(k,kp,k2,epsPsi);
                    multTable[3][3]=spp(k,kp,k2,k2);
                    multTable[3][4]=spp(k,kp,k2,k3);
                    multTable[4][0]=spp(k,kp,k3,epsG2);
                    multTable[4][1]=spp(k,kp,k3,epsG3);
                    multTable[4][2]=spp(k,kp,k3,epsPsi);
                    multTable[4][3]=spp(k,kp,k3,k2);
                    multTable[4][4]=spp(k,kp,k3,k3);
                    for(int iAmp=0; iAmp<nAmps; ++iAmp) {
                        for(int ciAmp=0; ciAmp<nAmps; ++ciAmp) {
                            matr2 += real(Amp[iAmp]*conj(Amp[ciAmp])*multTable[iAmp][ciAmp]);
                        };
                    };                    
                }
            }
        }
//        for(int _i=0; _i<4; ++_i) epsG2[_i]=k2[_i];
        
         
         cout<<" matr2="<<matr2<<endl;

        if(W2<minW2||W2>maxW2) continue;
        if(z<zMin || z>zMax) continue;
        if(Q2<minQ2||Q2>maxQ2) continue;
            
        Q2_scale=pow(xi*pT(pPsi),2);
        Q2_scale = xi*xi*(Q2+Mcc*Mcc);
        alphas=lhapdf_pdf->alphasQ2((double) Q2_scale);
        dbl_type pdf = lhapdf_pdf->xfxQ2(0, (double)x, (double)Q2_scale)/x;

        dbl_type wt=ramEP->wt;
        if(!(matr2>0)) {
            nNegative++;
            continue;
        };
        
        
        

        ++nPassed;
    };
 }