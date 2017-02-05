#include "kinematics/RamboEP.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include "kinematics/Rambo3.h"

dbl_type Pi = acos(-1), alpha = 1. / 137, alphaS = 0.3, ec = 2. / 3;
dbl_type M = 3.1;

dbl_type F(dbl_type s, dbl_type t, dbl_type u) {
    return (64 * alpha * pow(alphaS, 2) * pow(ec, 2) * pow(M, -1) * pow(Pi, 2) * pow(s, -4) * pow(s + t, -2) * pow(s + u, -2) * pow(t + u, -2)) / 27.;
}

dbl_type T(dbl_type Q, dbl_type s, dbl_type t, dbl_type u) {
    return -4 * pow(Q, 6) * pow(t, 2)*(pow(s, 2) + pow(t, 2)) + 2 * t * pow(Q, 4)*(3 * t * (t + u) * pow(s, 2) + (3 * t - 2 * u) * pow(s, 3) + 2 * s * (t - u) * pow(t, 2) + 2 * (t + u) * pow(t, 3)) -
            2 * s * pow(Q, 2)*(-2 * t * u * (t + u) * pow(s, 2) - s * (2 * t - u) * u * pow(t, 2) - 2 * u * (t + u) * pow(t, 3) + pow(s, 3) * pow(t - u, 2)) +
            2 * pow(s, 2)*((t + u) * pow(t, 2) * pow(u, 2) + pow(s, 3)*(t * u + pow(t, 2) + pow(u, 2)) + s * t * u * (3 * t * u + pow(t, 2) + pow(u, 2)) + pow(s, 2) * pow(t + u, 3));
}

dbl_type L(dbl_type Q, dbl_type s, dbl_type t, dbl_type u) {
    return 2 * t * (s * (t - u) + t * (t + u)) * pow(Q, 2)*(pow(s, 2) - 2 * pow(t, 2)) - 2 * pow(Q, 4)*(pow(s, 2) - 2 * pow(t, 2)) * pow(t, 2) -
            s * (2 * (t + u) * pow(s, 2) * pow(t, 2) + 4 * u * (t + u) * pow(t, 3) + pow(s, 3)*(pow(t, 2) + pow(u, 2)) + s * pow(t, 2)*(6 * t * u + pow(t, 2) + pow(u, 2)));
}

int main(void) {
    Random *random = new Random();
    dbl_type ecm = 600, S = pow(ecm, 2), x=1;
    Rambo3 ram3(0, 0, M, random);
    dbl_type P[4], p[4], k[4], kp[4], pPsi[4], q[4], kg[4];
    set_v4(P, 0, 0, -ecm / 2, ecm / 2);
    set_v4(p, 0, 0, -ecm / 2, ecm / 2);
    set_v4(k, 0, 0, ecm / 2, ecm / 2);
        assert(is_zero(mass2(P)));
//        println_4v("P", P);
        assert(is_zero(mass2(k)));
//        println_4v("k", k);
        assert(are_equal(sum_mass2(P, k), S));
        for(int iEv=0; iEv<1e6; ++iEv) {
            ram3.next(ecm, kp, kg, pPsi);
        //        println_4v("kp",kp);
        //        println_4v("kg",kg);
        //        println_4v("pPsi",pPsi);
                assert(is_zero(mass2(kp)));
                assert(is_zero(mass2(kg)));
                assert(are_equal(mass2(pPsi),M*M));
                dbl_type PtotOut[4];
                sum(kp,kg,pPsi,PtotOut);
        //        println_4v("PtotOut",PtotOut);
                assert(is_zero(PtotOut[0]));
                assert(is_zero(PtotOut[1]));
                assert(is_zero(PtotOut[2]));
                assert(are_equal(PtotOut[3],ecm));
            subtract(k,kp,q);
            dbl_type Q2=-mass2(q);
            dbl_type y=sp(q,P)/sp(k,P);
            dbl_type s=sum_mass2(q,p);
            dbl_type t=subtract_mass2(q,pPsi);
            dbl_type u=subtract_mass2(p,pPsi);
                assert(are_equal(s,x*y*S-Q2));
                assert(are_equal(s+t+u,M*M-Q2));
            dbl_type FF=F(s,t,u);
            dbl_type TT=T(sqrt(Q2),s,t,u);
            dbl_type LL=L(sqrt(Q2),s,t,u);
            dbl_type d3Sigma = -alpha/(2*Pi)*FF*((1+pow(1-y,2))/(y*Q2)*TT-4*(1-y)/y*LL);
//            cout<<" d3Sigma="<<d3Sigma<<endl;
            if(d3Sigma<0) cout<<" negative!!!"<<endl;
    };
    //    RamboEP ramEP(ecm, random, M, 0, 0);
    //    dbl_type P[4], k[4], kp[4], pPsi[4], q[4], kg[4];
    //    dbl_type minY = M * M / S, maxY = 1;
    //    ramEP.set_minmaxY(minY / S, maxY);
    //    set_v4(P,ramEP.P);
    //        assert(is_zero(mass2(P)));
    //        println_4v_math("P",P);
    //    set_v4(k,ramEP.kIn);
    //        assert(is_zero(mass2(k)));
    //        assert(are_equal(sum_mass2(P,k),S));
    //        println_4v_math("k",k);
    //        cout<<"-------"<<endl;
    //    dbl_type x=1;
    //    while(ramEP.next(kp, pPsi, kg, x)) {};
    //        println_4v_math("kp",kp);
    //        println_4v_math("pPsi",pPsi);
    //        println_4v_math("kg",kg);
    //        dbl_type Pg[4];
    //        set_v4(Pg,P);
    //        mult(Pg,1-x);
    //        println_4v_math("Pg",Pg);
    //        cout<<"-------"<<endl;
    //        assert(is_zero(mass2(kp)));
    //        assert(are_equal(mass2(pPsi),M*M));
    //        assert(is_zero(mass2(kg)));
    //        dbl_type PtotIn[4], PtotOut[4];
    //        sum(k,P,PtotIn);
    //        println_4v_math("PtotIn",PtotIn);
    //        sum(kp,pPsi,kg,PtotOut);
    //        println_4v_math("PtotOut",PtotOut);
    //
    //        assert(are_equal(PtotIn[0],PtotOut[0]));
    //        assert(are_equal(PtotIn[1],PtotOut[1]));
    //        assert(are_equal(PtotIn[2],PtotOut[2]));
    //        assert(are_equal(PtotIn[3],PtotOut[3]));
    return 0;
}

