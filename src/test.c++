#define CATCH_CONFIG_MAIN 
#include "catch.hpp"

#include "kinematics/algebra.h"
#include "kinematics/Rambo3.h"
#include "kinematics/Random.h"
#include "kinematics/Rambo2.h"
#include "kinematics/RamboEP.h"

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
    int nEv = 1e5;
    dbl_type sum = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        dbl_type wt = rambo.next(mmu, k1, k2, k3);
        dbl_type matr2 = 128 * sp(P, k1) * sp(k2, k3);
        sum += wt*matr2;
    };
    sum /= nEv;
    sum /= 2 * 2 * mmu;
    dbl_type th = pow(mmu, 5) / (192 * pow(PI, 3));
    //    cout<<" sum="<<sum<<" th="<<th<<" sum/th="<<sum/th<<endl;
    REQUIRE(sum == Approx(th));
}

TEST_CASE("ramboEP") {
    Random *random_generator = new Random();
    dbl_type ecm = 1e3, m1=10, m2=20;
    dbl_type P[4], k[4], kp[4], q[4], k1[4], k2[4], Ptot[4];

    RamboEP rambo(ecm, random_generator,m1, m2);
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
    REQUIRE( ecm*ecm+subtract_mass2(k,kp)+subtract_mass2(P,kp) == Approx(rambo.W2));
    REQUIRE(mass2(k1)==Approx(m1*m1));
    REQUIRE(mass2(k2)==Approx(m2*m2));
    
    sum(k1,k2,Ptot);
    REQUIRE(mass2(Ptot)==Approx(rambo.W2));
    sum(kp,k1,k2,Ptot);
    REQUIRE(mass2(Ptot)==Approx(ecm*ecm));
    REQUIRE(Ptot[0]==Approx(0));
    REQUIRE(Ptot[1]==Approx(0));
    REQUIRE(Ptot[2]==Approx(0));
    REQUIRE(Ptot[3]==Approx(ecm));    
}

TEST_CASE("EPflat") {
    dbl_type ecm=10;
    Random *random_generator = new Random();
    int nEv=1e5;
    dbl_type kP[4], k1[4], k2[4];
    dbl_type m1=random_generator->rand(0,ecm/2), m2=random_generator->rand(0,ecm/2);
    // 3D integration
    Rambo3 ram3(0,m1,m2,random_generator);
    dbl_type sum=0;
    for(int iEv=0; iEv<nEv; ++iEv) {
        sum += ram3.next(ecm,kP,k1,k2);
    };
    sum = sum/nEv;
    // eP integration
    RamboEP ramEP(ecm, random_generator,m1,m2);
    Rambo2 ram2(m1,m2,random_generator);
    dbl_type sumEP=0;
    for(int iEv=0; iEv<nEv; ++iEv) {
        ramEP.next(kP,k1,k2);
        sumEP += ramEP.wt;
    };
    sumEP = sumEP/nEv;
    REQUIRE(sum==Approx(sumEP).epsilon(1e-4));
}

TEST_CASE("toy1EP") {
    dbl_type ecm=100;
    Random *random_generator=new Random();
    dbl_type m1=random_generator->rand(0,ecm/2), m2=random_generator->rand(0,ecm/2);
//    dbl_type m1=ecm/3, m2=ecm/3;
    int nEv=1e4;
    dbl_type P[4], k[4], kp[4], k1[4], k2[4];
    set_v4(P,0,0,ecm/2,ecm/2);
    set_v4(k,0,0,-ecm/2,ecm/2);
    dbl_type s=ecm*ecm, t, u;
    dbl_type PI=acos(-1), alpha=1./137;
    // ram3
    Rambo3 ram3(0,m1,m2,random_generator);
    dbl_type sum=0;
    for(int iEv=0; iEv<nEv; ++iEv) {
        dbl_type wt=ram3.next(ecm,kp,k1,k2);
        u=subtract_mass2(k,kp);
        t=subtract_mass2(P,kp);
        dbl_type matr2 = -128*PI*PI*alpha*alpha*s*u/pow(s+u,2);
        sum += matr2*wt/nEv;
    }
//    cout<<" sum="<<sum<<endl;
    // ep
    RamboEP ramEP(ecm, random_generator, m1, m2);
    dbl_type sumEP=0;
    for(int iEv=0; iEv<nEv; ++iEv) {
        ramEP.next(kp,k1,k2);
        dbl_type matr2T=0;
        dbl_type matr2L = -4*PI*alpha*ramEP.Q2;
        sumEP += (matr2L*ramEP.wL+matr2T*ramEP.wT)/nEv;
    };
//    cout<<" sumEP="<<sumEP<<endl;
//    cout<<" sumEP/sum="<<sumEP/sum<<endl;
    REQUIRE(sumEP==Approx(sum).epsilon(0.001));
}
