#include "catch.hpp"

#include "kinematics/algebra.h"
#include "kinematics/Random.h"
#include "kinematics/Rambo2.h"

TEST_CASE("polarizations") {
    dbl_type M=3.1, P[4], k1[4], k2[4];
    cmplx eps[4];
    Random *random=new Random();
    set_v4(P,0,0,0,M);
    Rambo2 rambo(0,0,random);
    REQUIRE(mass2(P)==Approx(M*M));

    set_psi_polarization(0,P,eps); assert(is_zero(sp(P,eps)));
    set_psi_polarization(1,P,eps); REQUIRE(abs(sp(eps,P))==Approx(0));
    set_psi_polarization(2,P,eps); REQUIRE(abs(sp(eps,P))==Approx(0));

    rambo.next(M,k1,k2);
    REQUIRE(mass2(k1)==Approx(0));
    REQUIRE(mass2(k2)==Approx(0));
    set_gluon_polarization(0, k1, eps); REQUIRE(abs(sp(eps,k1))==Approx(0));
    set_gluon_polarization(1, k1, eps); REQUIRE(abs(sp(eps,k1))==Approx(0));
    set_gluon_polarization(0, k2, eps); REQUIRE(abs(sp(eps,k2))==Approx(0));
    set_gluon_polarization(1, k2, eps); REQUIRE(abs(sp(eps,k2))==Approx(0));
}

TEST_CASE("pol_sum") {
    dbl_type M=3.1, P[4], k1[4], k2[4];
    cmplx eps[4], ceps[4];
    Random *random=new Random();
    set_v4(P,0,0,0,M);
    Rambo2 rambo(0,0,random);
    REQUIRE(mass2(P)==Approx(M*M));
    cmplx sum=0;
    for(int i=0; i<3; ++i) {
        set_psi_polarization(i,P,eps, ceps);
        sum += sp(eps, ceps);
    };
    REQUIRE(sum.real()==Approx(-3));
    REQUIRE(sum.imag()==Approx(0));


    rambo.next(M,k1,k2);
    REQUIRE(mass2(k1)==Approx(0));
    REQUIRE(mass2(k2)==Approx(0));
    cmplx sum1=0, sum2=0;
    for(int i=0; i<2; ++i) {
        set_gluon_polarization(i,k1,eps,ceps);
        sum1 += sp(eps, ceps);
        set_gluon_polarization(i,k2,eps,ceps);
        sum2 += sp(eps, ceps);    
    };
    REQUIRE(sum1.real()==Approx(-2));
    REQUIRE(sum1.imag()==Approx(0));
    REQUIRE(sum2.real()==Approx(-2));
    REQUIRE(sum2.imag()==Approx(0));

}

TEST_CASE("lepton_current") {
    Random random(10);
    dbl_type k1[4], k2[4],q[4];
    cmplx eps[4];
    dbl_type w1=random.rand(0,1), w2=random.rand(0,1);
    dbl_type cos1=random.rand(-1,1), sin1=sqrt(1-cos1*cos1);
    set_v4(k1,0,w1,0,w1);
    REQUIRE(is_zero(mass2(k1)));
    dbl_type cos2=random.rand(-1,1), sin2=sqrt(1-cos2*cos2);
    set_v4(k2,0,w2*sin2, w2*cos2,w2);
    REQUIRE(is_zero(mass2(k2)));
    subtract(k1,k2,q);
    lepton_current(0,k1,k2,eps);
    REQUIRE(sp(q,eps).real()==Approx(0));
    //
    REQUIRE_FALSE(sp(eps,q).real()==Approx(1));
    REQUIRE(sp(eps,q).imag()==Approx(0));

    lepton_current(1,k1,k2,eps);
    REQUIRE(sp(eps,q).real()==Approx(0));
    REQUIRE_FALSE(sp(eps,q).real()==Approx(1));
    REQUIRE(sp(eps,q).imag()==Approx(0));
    
    // normalization
    cmplx ceps[4];
    lepton_current(0,k1,k2,eps);    conj_v4(ceps,eps);
    REQUIRE( sp(eps,ceps).real() == Approx(-4*sp(k1,k2)));
    REQUIRE( sp(eps,ceps).imag() == Approx(0));
    lepton_current(1,k1,k2,eps);    conj_v4(ceps,eps);
    REQUIRE( sp(eps,ceps).real() == Approx(-4*sp(k1,k2)));
    REQUIRE( sp(eps,ceps).imag() == Approx(0));
    
    
}
