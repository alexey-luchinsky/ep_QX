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
    int nEv=1e8, nPassed=0, nNegative=0;
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