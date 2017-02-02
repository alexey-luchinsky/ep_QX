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
dbl_type ecm=100, x;

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


// k, kp, k1, k2, k3, p
dbl_type getMatr2(dbl_type(&k)[4], dbl_type(&kp)[4], dbl_type(&k1)[4], dbl_type(&k2)[4], 
                dbl_type(&k3)[4], dbl_type(&p)[4]) {
    return (-524288*Opsi*pow(alpha,2)*pow(alphas,2)*pow(mc,-1)*pow(PI,4)*pow(sp(k1,k1),-2)*pow(sp(k2,p),-2)*pow(2*sp(k2,k3) + sp(k2,p) - sp(k3,p),-2)*pow(sp(k3,p),-2)*
     (-(pow(sp(k2,p),4)*pow(sp(k3,kp),2)) - pow(sp(k2,kp),2)*pow(sp(k3,p),3)*(4*pow(mc,2) + sp(k3,p)) + 
       pow(sp(k2,p),2)*sp(k3,p)*((pow(sp(k2,kp),2) - pow(sp(kp,p),2))*sp(k3,p) + 6*sp(k2,kp)*sp(k3,kp)*sp(k3,p) + pow(sp(k3,kp),2)*(-4*pow(mc,2) + sp(k3,p))) - 
       2*pow(sp(k3,p),2)*sp(k2,kp)*sp(k2,p)*(sp(k2,kp)*(-2*pow(mc,2) + sp(k3,p)) + sp(k3,p)*(sp(k3,kp) - sp(kp,p))) + 
       2*pow(sp(k2,p),3)*sp(k3,kp)*(sp(k3,kp)*(2*pow(mc,2) - sp(k3,p)) - sp(k3,p)*(sp(k2,kp) + sp(kp,p))) - 
       4*pow(mc,2)*pow(sp(k2,k3),3)*(-(sp(k3,kp)*sp(kp,p)) + sp(k2,kp)*(2*sp(k3,kp) + sp(kp,p))) + 
       pow(sp(k2,k3),2)*(pow(sp(k3,kp),2)*(32*pow(mc,4) - pow(sp(k2,p),2) - 4*pow(mc,2)*sp(k2,p)) + 
          pow(sp(k2,kp),2)*(32*pow(mc,4) - pow(sp(k3,p),2) + 4*pow(mc,2)*sp(k3,p)) + 2*pow(sp(kp,p),2)*(-4*pow(mc,4) + sp(k2,p)*sp(k3,p)) - 
          2*sp(k3,kp)*(sp(k2,p)*(4*pow(mc,2) - sp(k3,p)) + 2*pow(mc,2)*(4*pow(mc,2) + sp(k3,p)))*sp(kp,p) + 
          2*sp(k2,kp)*(2*sp(k3,kp)*(8*pow(mc,4) + sp(k2,p)*(pow(mc,2) - sp(k3,p)) - pow(mc,2)*sp(k3,p)) + 
             (8*pow(mc,4) - 4*pow(mc,2)*sp(k3,p) - sp(k2,p)*(2*pow(mc,2) + sp(k3,p)))*sp(kp,p))) + 
       2*sp(k2,k3)*(pow(sp(k2,kp),2)*sp(k3,p)*(sp(k2,p)*(-6*pow(mc,2) + sp(k3,p)) + sp(k3,p)*(4*pow(mc,2) + sp(k3,p))) + 
          sp(k2,kp)*(pow(sp(k2,p),2)*(sp(k3,kp)*(2*pow(mc,2) - 3*sp(k3,p)) - sp(k3,p)*sp(kp,p)) + 
             2*pow(mc,2)*sp(k3,p)*(sp(k3,kp)*(-4*pow(mc,2) + sp(k3,p)) + sp(k3,p)*sp(kp,p)) + 
             sp(k2,p)*(sp(k3,kp)*(8*pow(mc,4) + 3*pow(sp(k3,p),2) - 8*pow(mc,2)*sp(k3,p)) - 6*pow(mc,2)*sp(k3,p)*sp(kp,p))) - 
          sp(k2,p)*(pow(sp(k2,p),2)*pow(sp(k3,kp),2) + sp(k2,p)*(-(pow(sp(kp,p),2)*sp(k3,p)) + pow(sp(k3,kp),2)*(-4*pow(mc,2) + sp(k3,p)) + 
                2*pow(mc,2)*sp(k3,kp)*sp(kp,p)) + sp(k3,p)*(6*pow(mc,2)*pow(sp(k3,kp),2) + pow(sp(kp,p),2)*(-2*pow(mc,2) + sp(k3,p)) + 
                sp(k3,kp)*(-6*pow(mc,2) + sp(k3,p))*sp(kp,p))))))/81.;
}



int main(void) {
    string pdf_set_name="ct10";
    lhapdf_pdf=LHAPDF::mkPDF(pdf_set_name,0);
    Random *random=new Random();

    // experimental conditions from Kinehl
    ecm=300;
    RamboEP *ramEP=new RamboEP(ecm, random, Mcc,0);
    dbl_type minQ2=2, maxQ2=80;
    ramEP->set_minmaxQ2(minQ2, maxQ2);
    dbl_type minW2=pow(40,2), maxW2=pow(180,2);
    
    
    TNtuple tup("tup","tup","Q2:Y:pTpsi:W2:matr2:wt:pdf");
    TH1F hPT2("hPT2","hPT2",50,0,30); hPT2.Sumw2();
    TFile file("out.root","RECREATE");

    dbl_type Q2_scale=pow(Mcc,2);
    dbl_type k[4], kp[4], k1[4], k2[4], k3[4], p[4];
    int nEv=1e8, nPassed=0, nNegative=0;
    dbl_type sum=0, dsigma, sigma=0;
    for(int iEv=0; iEv<nEv; ++iEv) {
        if( iEv % (nEv/10) == 0) cout<<"---- Event "<<iEv<<" ("
                <<(int)(100.*iEv/nEv)<<" %) --- sigma="<<sigma*nEv/iEv<<" pb"<<endl;
        x=random->rand(0,1);
        if(!kinematics(ramEP, k, kp, k1, k2, k3, p)) continue;
        Q2_scale=2*pT_squared(p);
        alphas=lhapdf_pdf->alphasQ2((double) Q2_scale);
        dbl_type pdf = lhapdf_pdf->xfxQ2(0, (double)x, (double)Q2_scale)/x;

        dbl_type wt=ramEP->wt;
        dbl_type matr2=getMatr2(k, kp, k1, k2, k3, p);
        if(!(matr2>0)) {
            nNegative++;
            continue;
        };
        dbl_type W2=ramEP->W2;
        if(W2<minW2||W2>maxW2) continue;

        sum += matr2*wt*pdf/nEv;
        dsigma = picob*matr2*wt*pdf/(64*pow(ecm,2)*x*ramEP->Y)/nEv;
        sigma += dsigma;

        tup.Fill(ramEP->Q2,ramEP->Y,pT(p), W2, matr2,wt/nEv, pdf);
        hPT2.Fill(pT_squared(p),dsigma);
        ++nPassed;
    };
    tup.Write(); hPT2.Write();
    file.Save();
    write_histogram_to_file(hPT2,"hPT2.hst");
    cout<<" sigma="<<sigma<<" nb"<<endl;
    cout<<nPassed<<" ("<<(int)(100.*nPassed/nEv)<<"%) events passed"<<endl;
    cout<<nNegative<<" ("<<(int)(100.*nNegative/nEv)<<"%) events with negative matr2"<<endl;
    cout<<ramEP->nFault<<" faults in ramEP"<<endl;
}