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
