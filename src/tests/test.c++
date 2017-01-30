#define CATCH_CONFIG_MAIN 
#include "catch.hpp"

#include "kinematics/algebra.h"
#include "kinematics/Rambo3.h"
#include "kinematics/Random.h"
#include "kinematics/Rambo2.h"
#include "kinematics/RamboEP.h"

dbl_type PI = acos(-1), alpha = 1. / 137;
Random *random_generator = new Random();
dbl_type M, m1, m2, m3;
dbl_type P[4], k[4], kp[4], q[4], k1[4], k2[4], k3[4], Ptot[4];
int nEv = 1e4;
dbl_type sum_, sumEP_;






