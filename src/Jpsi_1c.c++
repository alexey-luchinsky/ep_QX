#include "kinematics/RamboEP.h"

int main(void) {
    dbl_type Mcc=3.1, mc=Mcc/2, ecm=10, x=0.1;
    Random *random=new Random();
    RamboEP ramEP(ecm, random, Mcc,0);
    dbl_type k[4], kp[4], k1[4], k2[4], k3[4], p[4];
    int nEv=100;
    for(int iEv=0; iEv<nEv; ++iEv) {
        while(!ramEP.next(kp,p,k3,x));
        set_v4(k,ramEP.kIn);                // k
            assert(is_zero(mass2(k)));       
        set_v4(kp, ramEP.kOut);             // kp
            assert(is_zero(mass2(kp)));
        set_v4(k2,ramEP.Pg);                // k2
            assert(are_equal(sum_mass2(p,k3),ramEP.W2));
        subtract(k,kp,k1);                  // k1=k-kp
            assert(are_equal(mass2(k1),-ramEP.Q2)); 
            assert(are_equal(mass2(p),Mcc*Mcc)); 
            assert(is_zero(mass2(k3)));
            assert(are_equal(sum_mass2(k,k2),x*ecm*ecm));

        dbl_type PI=acos(-1), alpha=1./137, alphas=0.3;
        dbl_type matr2=(8388608*pow(alpha,2)*pow(alphas,2)*pow(PI,4)*pow(sp(k1,k1),-2)*pow(sp(k2,p),-2)*pow(2*sp(k2,k3) + sp(k2,p) - sp(k3,p),-2)*pow(sp(k3,p),-2)*
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
        cout<<" matr2="<<matr2<<endl;
        assert(matr2>0);
    };
}