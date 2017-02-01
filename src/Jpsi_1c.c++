#include "kinematics/RamboEP.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
dbl_type Mcc=3.1, mc=Mcc/2, ecm=10, x=0.1;
dbl_type PI=acos(-1), alpha=1./137, alphas=0.3;

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
    return (8388608*pow(alpha,2)*pow(alphas,2)*pow(PI,4)*pow(sp(k1,k1),-2)*pow(sp(k2,p),-2)*pow(2*sp(k2,k3) + sp(k2,p) - sp(k3,p),-2)*pow(sp(k3,p),-2)*
         (pow(sp(k2,p),4)*pow(sp(k3,kp),2) + pow(sp(k2,kp),2)*pow(sp(k3,p),3)*(4*pow(mc,2) + sp(k3,p)) - 
           pow(sp(k2,p),2)*sp(k3,p)*((pow(sp(k2,kp),2) - pow(sp(kp,p),2))*sp(k3,p) + 6*sp(k2,kp)*sp(k3,kp)*sp(k3,p) + pow(sp(k3,kp),2)*(-4*pow(mc,2) + sp(k3,p))) + 
           2*pow(sp(k3,p),2)*sp(k2,kp)*sp(k2,p)*(sp(k2,kp)*(-2*pow(mc,2) + sp(k3,p)) + sp(k3,p)*(sp(k3,kp) - sp(kp,p))) + 
           2*pow(sp(k2,p),3)*sp(k3,kp)*(sp(k3,kp)*(-2*pow(mc,2) + sp(k3,p)) + sp(k3,p)*(sp(k2,kp) + sp(kp,p))) + 
           4*pow(mc,2)*pow(sp(k2,k3),3)*(-(sp(k3,kp)*sp(kp,p)) + sp(k2,kp)*(2*sp(k3,kp) + sp(kp,p))) + 
           pow(sp(k2,k3),2)*(pow(sp(k3,kp),2)*(-32*pow(mc,4) + pow(sp(k2,p),2) + 4*pow(mc,2)*sp(k2,p)) + 
              pow(sp(k2,kp),2)*(-32*pow(mc,4) + pow(sp(k3,p),2) - 4*pow(mc,2)*sp(k3,p)) + 2*pow(sp(kp,p),2)*(4*pow(mc,4) - sp(k2,p)*sp(k3,p)) + 
              2*sp(k3,kp)*(sp(k2,p)*(4*pow(mc,2) - sp(k3,p)) + 2*pow(mc,2)*(4*pow(mc,2) + sp(k3,p)))*sp(kp,p) + 
              sp(k2,kp)*(-4*sp(k3,kp)*(8*pow(mc,4) + sp(k2,p)*(pow(mc,2) - sp(k3,p)) - pow(mc,2)*sp(k3,p)) + 
                 2*(-8*pow(mc,4) + 4*pow(mc,2)*sp(k3,p) + sp(k2,p)*(2*pow(mc,2) + sp(k3,p)))*sp(kp,p))) - 
           2*sp(k2,k3)*(pow(sp(k2,kp),2)*sp(k3,p)*(sp(k2,p)*(-6*pow(mc,2) + sp(k3,p)) + sp(k3,p)*(4*pow(mc,2) + sp(k3,p))) + 
              sp(k2,kp)*(pow(sp(k2,p),2)*(sp(k3,kp)*(2*pow(mc,2) - 3*sp(k3,p)) - sp(k3,p)*sp(kp,p)) + 
                 2*pow(mc,2)*sp(k3,p)*(sp(k3,kp)*(-4*pow(mc,2) + sp(k3,p)) + sp(k3,p)*sp(kp,p)) + 
                 sp(k2,p)*(sp(k3,kp)*(8*pow(mc,4) + 3*pow(sp(k3,p),2) - 8*pow(mc,2)*sp(k3,p)) - 6*pow(mc,2)*sp(k3,p)*sp(kp,p))) - 
              sp(k2,p)*(pow(sp(k2,p),2)*pow(sp(k3,kp),2) + sp(k2,p)*(-(pow(sp(kp,p),2)*sp(k3,p)) + pow(sp(k3,kp),2)*(-4*pow(mc,2) + sp(k3,p)) + 
                    2*pow(mc,2)*sp(k3,kp)*sp(kp,p)) + sp(k3,p)*(6*pow(mc,2)*pow(sp(k3,kp),2) + pow(sp(kp,p),2)*(-2*pow(mc,2) + sp(k3,p)) + 
                    sp(k3,kp)*(-6*pow(mc,2) + sp(k3,p))*sp(kp,p))))))/9.;
}



int main(void) {
    TNtuple tup("tup","tup","Q2:Y:pTpsi:matr2:wt");
    Random *random=new Random();
    RamboEP *ramEP=new RamboEP(ecm, random, Mcc,0);
    ramEP->minQ2=0.3; ramEP->maxQ2=0.5;
    string dir="Q2_"+std::to_string(ramEP->minQ2)+"_"+std::to_string(ramEP->maxQ2);
    system(("mkdir "+dir).c_str());
    TFile file((dir+"/out.root").c_str(),"RECREATE");
    TH1F hQ2("hQ2","hQ2",20,ramEP->minQ2,ramEP->maxQ2);
    TH1F hQ2M("hQ2M","hQ2M",20,ramEP->minQ2,ramEP->maxQ2);
    dbl_type k[4], kp[4], k1[4], k2[4], k3[4], p[4];
    int nEv=1e7, nPassed=0, nNegative=0;
    dbl_type sum=0;
    for(int iEv=0; iEv<nEv; ++iEv) {
        if( iEv % (nEv/10) == 0) cout<<"---- Event "<<iEv<<" ("<<(int)(100.*iEv/nEv)<<" %) -----"<<endl;
        kinematics(ramEP, k, kp, k1, k2, k3, p);

        dbl_type wt=ramEP->wt;
        dbl_type matr2=getMatr2(k, kp, k1, k2, k3, p);
        if(!(matr2>0)) {
//            cout<<"======= matr2<0 ======="<<endl;
//            println_4v("k",k);
//            println_4v("kp",kp);
//            println_4v("k1",k1);
//            println_4v("k2",k2);
//            println_4v("k3",k3);
//            println_4v("p",p);
//            cout<<" W2="<<ramEP->W2<<endl;
//            cout<<" Q2="<<ramEP->Q2<<endl;
//            cout<<" Y="<<ramEP->Y<<endl;
//            cout<<" wt="<<wt<<endl;
//            cout<<"matr2:"<<matr2<<endl;
            nNegative++;
            continue;
        };
//        assert(matr2>0);
        tup.Fill(ramEP->Q2,ramEP->Y,pT(p), matr2,wt/nEv);
        hQ2M.Fill(ramEP->Q2,matr2*wt/nEv);
        hQ2.Fill(ramEP->Q2,wt/nEv);
        sum += matr2*wt/nEv;
        ++nPassed;
    };
    tup.Write(); hQ2.Write(); hQ2M.Write();
    file.Save();
    write_histogram_to_file(hQ2, (dir+"/hQ2.txt").c_str());
    write_histogram_to_file(hQ2M, (dir+"/hQ2M.txt").c_str());

    cout<<" minQ2="<<ramEP->minQ2<<" maxQ2="<<ramEP->maxQ2<<endl;
    cout<<" sum="<<sum<<endl;
    cout<<nPassed<<" ("<<(int)(100.*nPassed/nEv)<<"%) events passed"<<endl;
    cout<<nNegative<<" ("<<(int)(100.*nNegative/nEv)<<"%) events with negative matr2"<<endl;
    cout<<ramEP->nFault<<" faults in ramEP"<<endl;
}