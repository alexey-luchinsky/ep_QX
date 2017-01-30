#include "catch.hpp"

#include "kinematics/algebra.h"
#include "kinematics/Rambo3.h"
#include "kinematics/Random.h"
#include "kinematics/Rambo2.h"
#include "kinematics/RamboEP.h"

extern dbl_type PI, alpha, M, m1, m2, m3, sum_, sumEP_;
extern int nEv;
extern Random *random_generator;
extern dbl_type P[4], k[4], kp[4], q[4], k1[4], k2[4], k3[4], Ptot[4];

TEST_CASE("toy1EP") {
    dbl_type ecm = 10;
    dbl_type m1 = random_generator->rand(0, ecm / 2), m2 = random_generator->rand(0, ecm / 2);
    int nEv = 1e6;
    set_v4(P, 0, 0, ecm / 2, ecm / 2);
    set_v4(k, 0, 0, -ecm / 2, ecm / 2);
    dbl_type s = ecm*ecm, t, u;
    dbl_type PI = acos(-1), alpha = 1. / 137;
    // ram3
    Rambo3 ram3(0, m1, m2, random_generator);
    dbl_type sum = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        dbl_type wt = ram3.next(ecm, kp, k1, k2);
        u = subtract_mass2(k, kp);
        t = subtract_mass2(P, kp);
        dbl_type matr2 = -128 * PI * PI * alpha * alpha * s * u / pow(s + u, 2);
        sum += matr2 * wt / nEv;
    }
    // ep
    RamboEP ramEP(ecm, random_generator, m1, m2);
    dbl_type sumEP = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        ramEP.next(kp, k1, k2);
        dbl_type matr2T = 0;
        dbl_type matr2L = -4 * PI * alpha * ramEP.Q2;
        sumEP += (matr2L * ramEP.wL + matr2T * ramEP.wT) / nEv;
    };
    REQUIRE(sumEP / sum == Approx(1).epsilon(0.01));
};

TEST_CASE("toy2EP") {
    dbl_type ecm = 3;
    dbl_type m1 = random_generator->rand(0, 0.1 * ecm), m2 = random_generator->rand(0, 0.1 * ecm);
    //    m1=0; m2=0;
    int nEv = 1e5;
    set_v4(P, 0, 0, -ecm / 2, ecm / 2);
    set_v4(k, 0, 0, ecm / 2, ecm / 2);
    dbl_type s = ecm*ecm, t, u;
    // ram3
    Rambo3 ram3(0, m1, m2, random_generator);
    dbl_type sum = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        dbl_type wt = ram3.next(ecm, kp, k1, k2);
        subtract(k, kp, q);
        t = subtract_mass2(k, kp);
        u = subtract_mass2(P, kp);
        dbl_type Q2 = -t;
        //        dbl_type matr2=1;
        dbl_type matr2 = 32 * pow(alpha, 2) * pow(PI, 2)*(pow(m1, 2) * pow(m2, 2) * pow(Q2, 2) + Q2 * pow(m2, 2) * pow(sp(k1, q), 2) + Q2 * pow(m1, 2) * pow(sp(k2, q), 2) - 4 * pow(m2, 2) * sp(k, q) * sp(k1, kp) * sp(k1, q) -
                4 * Q2 * pow(m1, 2) * sp(k, k2) * sp(k2, kp) - 4 * pow(sp(k1, q), 2) * sp(k, k2) * sp(k2, kp) + 4 * sp(k, k2) * sp(k1, kp) * sp(k1, q) * sp(k2, q) - 4 * pow(m1, 2) * sp(k, q) * sp(k2, kp) * sp(k2, q) +
                4 * pow(m1, 2) * pow(m2, 2) * sp(k, q) * sp(kp, q) - 4 * pow(m1, 2) * sp(k, k2) * sp(k2, q) * sp(kp, q) - pow(sp(k1, k2), 2)*(pow(Q2, 2) + 4 * sp(k, q) * sp(kp, q)) -
                4 * sp(k, k1)*((Q2 * pow(m2, 2) + pow(sp(k2, q), 2)) * sp(k1, kp) + sp(k1, q)*(-(sp(k2, kp) * sp(k2, q)) + pow(m2, 2) * sp(kp, q))) +
                sp(k1, k2)*(4 * sp(k, q) * sp(k1, q) * sp(k2, kp) + 4 * sp(k, q) * sp(k1, kp) * sp(k2, q) - 2 * Q2 * sp(k1, q) * sp(k2, q) + 4 * sp(k, k2)*(Q2 * sp(k1, kp) + sp(k1, q) * sp(kp, q)) +
                4 * sp(k, k1)*(Q2 * sp(k2, kp) + sp(k2, q) * sp(kp, q))));

        sum += matr2 * wt / nEv;
    }
    //    cout<<" sum="<<sum<<endl;

    // ep
    RamboEP ramEP(ecm, random_generator, m1, m2);
    dbl_type p[4];
    set_v4(p, ramEP.P);
    dbl_type sumEP = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        ramEP.next(kp, k1, k2);
        set_v4(q, ramEP.q);
        dbl_type Q2 = ramEP.Q2, W2 = ramEP.W2;

        dbl_type matr2T = -4 * alpha * PI * pow(Q2, 2) * pow(sp(p, q), -2)*(Q2 * (Q2 * pow(m1, 2) + pow(sp(k1, q), 2)) * pow(sp(k2, p), 2) + Q2 * pow(sp(k1, p), 2)*(Q2 * pow(m2, 2) + pow(sp(k2, q), 2)) +
                pow(sp(p, q), 2)*(pow(m2, 2) * pow(sp(k1, q), 2) + pow(m1, 2) * pow(sp(k2, q), 2) - 2 * sp(k1, k2) * sp(k1, q) * sp(k2, q)) +
                2 * Q2 * sp(k2, p)*(-(sp(k1, k2) * sp(k1, q)) + pow(m1, 2) * sp(k2, q)) * sp(p, q) -
                2 * Q2 * sp(k1, p)*(sp(k1, q)*(sp(k2, p) * sp(k2, q) - pow(m2, 2) * sp(p, q)) + sp(k1, k2)*(Q2 * sp(k2, p) + sp(k2, q) * sp(p, q))));

        dbl_type matr2L = 4 * alpha * PI * pow(Q2, 3) * pow(sp(p, q), -2)*((Q2 * pow(m1, 2) + pow(sp(k1, q), 2)) * pow(sp(k2, p), 2) + pow(sp(k1, p), 2)*(Q2 * pow(m2, 2) + pow(sp(k2, q), 2)) +
                (-(pow(m1, 2) * pow(m2, 2)) + pow(sp(k1, k2), 2)) * pow(sp(p, q), 2) + 2 * sp(k2, p)*(-(sp(k1, k2) * sp(k1, q)) + pow(m1, 2) * sp(k2, q)) * sp(p, q) -
                2 * sp(k1, p)*(sp(k1, q)*(sp(k2, p) * sp(k2, q) - pow(m2, 2) * sp(p, q)) + sp(k1, k2)*(Q2 * sp(k2, p) + sp(k2, q) * sp(p, q))));

        sumEP += (matr2L * ramEP.wL + matr2T * ramEP.wT) / nEv;
    };
    REQUIRE(sumEP / sum == Approx(1).epsilon(0.01));
}

TEST_CASE("Q2_cut") {
    dbl_type ecm = 100;
    dbl_type minQ2 = pow(ecm, 2) / 3, maxQ2 = 2 * pow(ecm, 2) / 3;
    dbl_type m1 = random_generator->rand(0, ecm / 2), m2 = random_generator->rand(0, ecm / 2);
    int nEv = 1e6;
    set_v4(P, 0, 0, ecm / 2, ecm / 2);
    set_v4(k, 0, 0, -ecm / 2, ecm / 2);
    dbl_type s = ecm*ecm, t, u, Q2;
    dbl_type PI = acos(-1), alpha = 1. / 137;
    // ram3
    Rambo3 ram3(0, m1, m2, random_generator);
    dbl_type sum = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        dbl_type wt = ram3.next(ecm, kp, k1, k2);
        u = subtract_mass2(k, kp);
        t = subtract_mass2(P, kp);
        Q2 = -t;
        if (Q2 < minQ2 || Q2 > maxQ2) continue;
        dbl_type matr2 = -128 * PI * PI * alpha * alpha * s * u / pow(s + u, 2);
        sum += matr2 * wt / nEv;
    }
    // ep
    RamboEP ramEP(ecm, random_generator, m1, m2);
    ramEP.minQ2 = minQ2;
    ramEP.maxQ2 = maxQ2;
    dbl_type sumEP = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        if (!ramEP.next(kp, k1, k2)) continue;
        ;
        dbl_type matr2T = 0;
        dbl_type matr2L = -4 * PI * alpha * ramEP.Q2;
        sumEP += (matr2L * ramEP.wL + matr2T * ramEP.wT) / nEv;
    };
    REQUIRE(sumEP / sum == Approx(1).epsilon(0.01));

}



TEST_CASE("toy1EPX") {
    dbl_type ecm = 10;
    dbl_type x = 0.5, ecm_ = sqrt(x) * ecm;
    dbl_type m1 = random_generator->rand(0, ecm_ / 2), m2 = random_generator->rand(0, ecm_ / 2);
    int nEv = 1e6;
    set_v4(k, 0, 0, ecm / 2, ecm / 2);
    set_v4(P, 0, 0, -ecm / 2, ecm / 2);

    dbl_type s = ecm*ecm, t, u;
    // ram3
    Rambo3 ram3(0, m1, m2, random_generator);
    dbl_type Pg[4];
    set_v4(Pg, 0, 0, -x * ecm / 2, x * ecm / 2);
    dbl_type Phard[4];
    sum(Pg, k, Phard);
    dbl_type sum = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        dbl_type wt = ram3.next(ecm_, kp, k1, k2);
        apply_boost_to(Phard, k1);
        apply_boost_to(Phard, k2);
        apply_boost_to(Phard, kp);

        dbl_type W2 = sum_mass2(k1, k2);
        dbl_type s = sum_mass2(k, Pg);
        dbl_type t = subtract_mass2(k, kp);
        dbl_type u = subtract_mass2(Pg, kp);

        dbl_type matr2 = -128 * PI * PI * alpha * alpha * s * u / pow(s + u, 2);
        sum += matr2 * wt / nEv;
    }
    // ep
    RamboEP ramEP(ecm, random_generator, m1, m2);
    dbl_type sumEP = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        ramEP.next(kp, k1, k2,x);
        dbl_type matr2T = 0;
        dbl_type matr2L = -4 * PI * alpha * ramEP.Q2;
        sumEP += (matr2L * ramEP.wL + matr2T * ramEP.wT) / nEv;
    };
    REQUIRE(sumEP / sum == Approx(1).epsilon(0.01));
};


