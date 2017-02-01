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

TEST_CASE("ramboEP") {
    dbl_type ecm = 1e3, m1 = 10, m2 = 20;

    RamboEP rambo(ecm, random_generator, m1, m2);
    set_v4(P, rambo.P);
    set_v4(k, rambo.kIn);
    REQUIRE(mass2(P) == Approx(0));
    REQUIRE(mass2(k) == Approx(0));
    REQUIRE(sp(P, k) == Approx(ecm * ecm / 2));

    rambo.next(kp, k1, k2);
    subtract(k, kp, q); // q = k-kp
    REQUIRE(mass2(kp) == Approx(0));
    REQUIRE(mass2(q) == Approx(-rambo.Q2));
    REQUIRE(sp(q, P) / sp(k, P) == Approx(rambo.Y));
    // s+t+u == W2
    REQUIRE(ecm * ecm + subtract_mass2(k, kp) + subtract_mass2(P, kp) == Approx(rambo.W2));
    REQUIRE(mass2(k1) == Approx(m1 * m1));
    REQUIRE(mass2(k2) == Approx(m2 * m2));

    sum(k1, k2, Ptot);
    REQUIRE(mass2(Ptot) == Approx(rambo.W2));
    sum(kp, k1, k2, Ptot);
    REQUIRE(mass2(Ptot) == Approx(ecm * ecm));
    REQUIRE(Ptot[0] == Approx(0));
    REQUIRE(Ptot[1] == Approx(0));
    REQUIRE(Ptot[2] == Approx(0));
    REQUIRE(Ptot[3] == Approx(ecm));
}

TEST_CASE("EPflat") {
    dbl_type ecm = 10;
    dbl_type m1 = random_generator->rand(0, ecm / 2), m2 = random_generator->rand(0, ecm / 2);
    nEv=1e5;
    // 3D integration
    Rambo3 ram3(0, m1, m2, random_generator);
    dbl_type sum = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        sum += ram3.next(ecm, kp, k1, k2);
    };
    sum = sum / nEv;
    // eP integration
    RamboEP ramEP(ecm, random_generator, m1, m2);
    Rambo2 ram2(m1, m2, random_generator);
    dbl_type sumEP = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        if(!ramEP.next(kp, k1, k2)) continue;
        sumEP += ramEP.wt;
    };
    sumEP = sumEP / nEv;
    REQUIRE(sumEP / sum == Approx(1).epsilon(0.01));
}

TEST_CASE("ramEPX_kin") {
    dbl_type ecm = 10, x = 0.5;
    RamboEP ramEP(ecm, random_generator);
    ramEP.next(kp, k1, k2, x);
    dbl_type Pg[4];
    set_v4(P, ramEP.P);
    set_v4(k, ramEP.kIn);
    set_v4(Pg, ramEP.Pg);
    set_v4(q, ramEP.q);
    REQUIRE(mass2(Pg) == Approx(0));
    REQUIRE(Pg[0] == Approx(x * P[0]));
    REQUIRE(Pg[1] == Approx(x * P[1]));
    REQUIRE(Pg[2] == Approx(x * P[2]));
    REQUIRE(Pg[3] == Approx(x * P[3]));
    REQUIRE(mass2(q) == Approx(-ramEP.Q2));
    REQUIRE(sp(P, q) / sp(P, k) == Approx(ramEP.Y));
    REQUIRE(sum_mass2(q, Pg) == Approx(ramEP.W2));
    dbl_type PhardIn[4], PhardOut[4];
    sum(Pg, q, PhardIn); // PhardIn = Pg + q
    sum(k1, k2, PhardOut); // PhardOut = k1+k2
    REQUIRE(mass2(PhardIn) == Approx(ramEP.W2));
    REQUIRE(mass2(PhardOut) == Approx(ramEP.W2));
    REQUIRE(PhardIn[0] == Approx(PhardOut[0]));
    REQUIRE(PhardIn[1] == Approx(PhardOut[1]));
    REQUIRE(PhardIn[2] == Approx(PhardOut[2]));
    REQUIRE(PhardIn[3] == Approx(PhardOut[3]));
}

TEST_CASE("EPflatX") {
    dbl_type ecm = 10;
    dbl_type x = 0.5;
    int nEv = 1e5;
    dbl_type m1 = random_generator->rand(0, ecm * sqrt(x) / 2), m2 = random_generator->rand(0, ecm * sqrt(x) / 2);
    // 3D integration
    Rambo3 ram3(0, m1, m2, random_generator);
    dbl_type sum = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        sum += ram3.next(sqrt(x) * ecm, kp, k1, k2);
    };
    sum = sum / nEv;
    // eP integration
    RamboEP ramEP(ecm, random_generator, m1, m2);
    Rambo2 ram2(m1, m2, random_generator);
    dbl_type sumEP = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        if(!ramEP.next(kp, k1, k2, x)) continue;
        sumEP += ramEP.wt;
    };
    sumEP = sumEP / nEv;
    REQUIRE(sum/sumEP == Approx(1).epsilon(1e-2));
}

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
    dbl_type ecm = 10;
    dbl_type x = 0.5, ecm_=sqrt(x)*ecm, s_=ecm_*ecm_;
    
    int nEv = 1e5;
    dbl_type m1 = random_generator->rand(0, ecm_/ 2), 
            m2 = random_generator->rand(0, ecm_ / 2);

    dbl_type minQ2=0.1*s_, midQ2=0.5*s_, maxQ2=0.9*s_;
    RamboEP ramEP_all(ecm, random_generator, m1, m2); ramEP_all.set_minmaxQ2(minQ2,maxQ2);
    RamboEP ramEP_left(ecm, random_generator, m1, m2); ramEP_left.set_minmaxQ2(minQ2, midQ2);
    RamboEP ramEP_right(ecm, random_generator, m1, m2); ramEP_right.set_minmaxQ2(midQ2, maxQ2);
    dbl_type sum_all = 0, sum_left=0, sum_right=0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        if(!ramEP_all.next(kp, k1, k2, x)) continue;
        if(!ramEP_left.next(kp, k1, k2, x)) continue;
        if(!ramEP_right.next(kp, k1, k2, x)) continue;
        sum_all += ramEP_all.wt/nEv;
        sum_left += ramEP_left.wt/nEv;
        sum_right += ramEP_right.wt/nEv;
    };
    dbl_type sum_sum=sum_left+sum_right;
    cout<<"[Q2cut_flat_sum] sum_all="<<sum_all<<" sum_sum="<<sum_sum<<" sum_all/sum_sum="<<sum_all/sum_sum<<endl;
    REQUIRE(sum_all/sum_sum==Approx(sum_all/sum_sum).epsilon(0.01));
}
