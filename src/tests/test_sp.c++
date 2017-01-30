#include "catch.hpp"


#include "kinematics/algebra.h"
#include "kinematics/Rambo3.h"
#include "kinematics/Random.h"

extern dbl_type PI, alpha, M, m1, m2, m3, sum_, sumEP_;
extern int nEv;
extern Random *random_generator;
extern dbl_type P[4], k[4], kp[4], q[4], k1[4], k2[4], k3[4], Ptot[4];

TEST_CASE("sp") {
    M = 3.1;
    set_v4(P, 0, 0, 0, M);
    set_v4(k, 0, 0, M, M);
    REQUIRE(sp(P, P) == Approx(pow(M, 2)));
    REQUIRE(sp(k, k) == Approx(0));
    REQUIRE(sp(P, k) == Approx(pow(M, 2)));
}

