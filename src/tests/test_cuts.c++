#include "catch.hpp"

#include "kinematics/algebra.h"
#include "kinematics/Rambo3.h"
#include "kinematics/Random.h"
#include "kinematics/Rambo2.h"
#include "kinematics/RamboEP.h"
#include "RamboEP.h"

extern dbl_type PI, alpha, M, m1, m2, m3, sum_, sumEP_;
extern int nEv;
extern Random *random_generator;
extern dbl_type P[4], k[4], kp[4], q[4], k1[4], k2[4], k3[4], Ptot[4];

TEST_CASE("Q2cut_flat") {
    dbl_type ecm = 10;
    dbl_type x = 0.5, ecm_=sqrt(x)*ecm;
    dbl_type minQ2 = pow(ecm_, 2) / 3, maxQ2 = 2 * pow(ecm_, 2) / 3;    
    int nEv = 1e5;
    dbl_type m1 = random_generator->rand(0, ecm_/ 2), 
            m2 = random_generator->rand(0, ecm_ / 2);
    set_v4(P, 0, 0, ecm_ / 2, ecm_ / 2);
    set_v4(k, 0, 0, -ecm_ / 2, ecm_ / 2);

    // 3D integration
    Rambo3 ram3(0, m1, m2, random_generator);
    dbl_type sum = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        dbl_type wt = ram3.next(ecm_, kp, k1, k2);
        dbl_type t = subtract_mass2(k, kp);
        dbl_type u = subtract_mass2(P, kp);
        dbl_type Q2 = -t;
        if (Q2 < minQ2 || Q2 > maxQ2) continue;
        sum += wt/nEv;
    };
    // eP integration
    RamboEP ramEP(ecm, random_generator, m1, m2);
    ramEP.set_minQ2(minQ2); ramEP.set_maxQ2(maxQ2);
    dbl_type sumEP = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        if(!ramEP.next(kp, k1, k2, x)) continue;
        sumEP += ramEP.wt/nEv;
    };
    cout<<"Q2cut_flat: sum="<<sum<<" sumEP="<<sumEP<<" sum/sumEP="<<sum/sumEP<<endl;
    REQUIRE(sum/sumEP == Approx(1).epsilon(1e-2));
}

TEST_CASE("Q2cut_flat_sum") {
    dbl_type ecm = 100;
    dbl_type x = 0.5, ecm_=sqrt(x)*ecm, s_=ecm_*ecm_;
    
    int nEv = 1e7;
    dbl_type m1 = random_generator->rand(0, ecm_/ 2), 
            m2 = random_generator->rand(0, ecm_ / 2);
    dbl_type minQ2=0, midQ2=0.5*s_, maxQ2=s_;
    RamboEP ramEP_all(ecm, random_generator, m1, m2); ramEP_all.set_minmaxQ2(minQ2,maxQ2);
    RamboEP ramEP_left(ecm, random_generator, m1, m2); ramEP_left.set_minmaxQ2(minQ2, midQ2);
    RamboEP ramEP_right(ecm, random_generator, m1, m2); ramEP_right.set_minmaxQ2(midQ2, maxQ2);
    dbl_type sum_all = 0, sum_left=0, sum_right=0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        ramEP_all.next(x);  ramEP_left.next(x); ramEP_right.next(x);
        sum_all += ramEP_all.wt/nEv;
        sum_left += ramEP_left.wt/nEv;
        sum_right += ramEP_right.wt/nEv;
    };
    dbl_type sum_sum=sum_left+sum_right;
    cout<<"[Q2cut_flat_sum] sum_all="<<sum_all<<" sum_sum="<<sum_sum<<" sum_all/sum_sum="<<sum_all/sum_sum<<endl;
    cout<<"[Q2cut_flat_sum] ramEP_all.nFault="<<ramEP_all.nFault<<
            " ramEP_left.nFault="<<ramEP_left.nFault<<
            " ramEP_right.nFault="<<ramEP_right.nFault<<endl;
    REQUIRE(sum_all/sum_sum==Approx(1).epsilon(0.01));
}

