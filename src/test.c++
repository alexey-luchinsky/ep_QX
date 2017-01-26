#define CATCH_CONFIG_MAIN 
#include "catch.hpp"

#include "kinematics/algebra.h"
#include "kinematics/Rambo3.h"
#include "kinematics/Random.h"
#include "kinematics/Rambo2.h"

dbl_type PI = acos(-1);

TEST_CASE("sp") {
    dbl_type M = 3.1;
    dbl_type P[4] = {0, 0, 0, M}, k[4] = {0, 0, M, M};
    REQUIRE(sp(P, P) == Approx(pow(M, 2)));
    REQUIRE(sp(k, k) == Approx(0));
    REQUIRE(sp(P, k) == Approx(pow(M, 2)));
}

TEST_CASE("rambo3") {
    Random *random_generator = new Random();
    dbl_type M = 1, m1 = 0.1, m2 = 0.2, m3 = 0.3;
    Rambo3 rambo(m1, m2, m3, random_generator);
    dbl_type k1[4], k2[4], k3[4], P[4];
    dbl_type wt = rambo.next(M, k1, k2, k3);
    REQUIRE(sp(k1, k1) == Approx(m1 * m1));
    REQUIRE(sp(k2, k2) == Approx(m2 * m2));
    REQUIRE(sp(k3, k3) == Approx(m3 * m3));
    sum(k1, k2, P);
    sum(k3, P, P);
    REQUIRE(sp(P, P) == Approx(M * M));
    REQUIRE(P[0] == Approx(0));
    REQUIRE(P[1] == Approx(0));
    REQUIRE(P[2] == Approx(0));
    REQUIRE(P[3] == Approx(M));
}

TEST_CASE("rambo2") {
    Random *random_generator = new Random();
    dbl_type M = 1, m1 = 0.1, m2 = 0.2;
    Rambo2 rambo(m1, m2, random_generator);
    dbl_type k1[4], k2[4], P[4];
    dbl_type WT = rambo.next(M, k1, k2);
    REQUIRE(sp(k1, k1) == Approx(m1 * m1));
    REQUIRE(sp(k2, k2) == Approx(m2 * m2));
    REQUIRE(sp(k1, k2) == Approx((M * M - m1 * m1 - m2 * m2) / 2));
    sum(k1, k2, P);
    REQUIRE(sp(P, P) == Approx(M * M));
    REQUIRE(P[0] == Approx(0));
    REQUIRE(P[1] == Approx(0));
    REQUIRE(P[2] == Approx(0));
    REQUIRE(P[3] == Approx(M));
    dbl_type lam = sqrt((1 - pow(m1 + m2, 2) / M / M)*(1 - pow(m1 - m2, 2) / M / M));
    REQUIRE(WT == Approx(lam / 8 / PI));
}

TEST_CASE("muon") {
    Random *random_generator = new Random();
    dbl_type mmu = 1;
    Rambo3 rambo(0, 0, 0, random_generator);
    dbl_type P[4] = {0, 0, 0, mmu}, k1[4], k2[4], k3[4];
    int nEv=1e5;
    dbl_type sum=0;
    for(int iEv=0; iEv<nEv; ++iEv) {
        dbl_type wt=rambo.next(mmu, k1, k2, k3);
        dbl_type matr2=128*sp(P,k1)*sp(k2,k3);
        sum += wt*matr2;
    };
    sum /= nEv;
    sum /= 2*2*mmu;
    dbl_type th=pow(mmu,5)/(192*pow(PI,3));
//    cout<<" sum="<<sum<<" th="<<th<<" sum/th="<<sum/th<<endl;
    REQUIRE( sum == Approx(th ));
    
}