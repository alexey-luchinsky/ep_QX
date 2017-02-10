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

TEST_CASE("eps") {
    dbl_type v[4], a[4], b[4], c[4], d[4];
    Random random;
    for(int i=0; i<4; ++i) {
        v[i]=random.rand(0,1);
        a[i]=random.rand(0,1);
        b[i]=random.rand(0,1);
        c[i]=random.rand(0,1);
        d[i]=random.rand(0,1);
    }
    REQUIRE(lcv(a,b,c,a)==Approx(0));
    REQUIRE(lcv(a,b,c,d)+lcv(b,a,c,d)==Approx(0));
    dbl_type tmp[4], tmp2[4], schouten[4];
    set_v4(schouten,v); mult(schouten,lcv(a,b,c,d));
    set_v4(tmp, a); mult(tmp, -lcv(v,b,c,d)); sum(schouten,tmp, tmp2); set_v4(schouten,tmp2);
    set_v4(tmp, b); mult(tmp, -lcv(a,v,c,d)); sum(schouten,tmp, tmp2); set_v4(schouten,tmp2);
    set_v4(tmp, c); mult(tmp, -lcv(a,b,v,d)); sum(schouten,tmp, tmp2); set_v4(schouten,tmp2);
    set_v4(tmp, d); mult(tmp, -lcv(a,b,c,v)); sum(schouten,tmp, tmp2); set_v4(schouten,tmp2);
    REQUIRE(is_zero(schouten[0]));
    REQUIRE(is_zero(schouten[1]));
    REQUIRE(is_zero(schouten[2]));
    REQUIRE(is_zero(schouten[3]));
}
