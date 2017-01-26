/* 
 * File:   rambo3.cpp
 * Author: luchinsky
 * 
 * Created on March 21, 2016, 11:34 AM
 */

#include "Rambo3.h"
#include "algebra.h"

// 3D LIPS q k3 -> k1 k2 k3

Rambo3::Rambo3(dbl_type m1_, dbl_type m2_, dbl_type m3_, Random *random_) {
    random = random_;
    m1 = m1_;
    m2 = m2_;
    m3 = m3_;
    min_cos = -1;
    max_cos = 1;
    PI = acos(-1);
}

dbl_type Rambo3::next(dbl_type ecm, dbl_type (&k1)[4], dbl_type (&k2)[4], dbl_type (&k3)[4]) {
    dbl_type min_q2 = pow(m1 + m2, 2);
    dbl_type max_q2 = pow(ecm - m3, 2);

    dbl_type WT = 1;

    // q k3 phase space
    dbl_type q2 = random->rand(min_q2, max_q2);
    WT *= (max_q2 - min_q2);
    dbl_type cosT = random->rand(min_cos, max_cos);
    dbl_type sinT = sqrt(1. - cosT * cosT);
    WT *= (max_cos - min_cos);

    dbl_type phi = random->rand(0, 2 * PI);
    dbl_type sinPhi = sin(phi);
    dbl_type cosPhi = cos(phi);

    WT *= 2 * PI;
    dbl_type eQ = (ecm * ecm + q2 - m3 * m3) / (2 * ecm);
    dbl_type e3 = (ecm * ecm + m3 * m3 - q2) / (2 * ecm);
    dbl_type mom = sqrt(eQ * eQ - q2);
    q[3] = eQ;
    q[0] = mom * sinT * cosPhi;
    q[1] = mom * sinT * sinPhi;
    q[2] = mom * cosT;

    k3[3] = e3;
    k3[0] = -mom * sinT * cosPhi;
    k3[1] = -mom * sinT * sinPhi;
    k3[2] = -mom * cosT;

    WT *= mom / (ecm * 4 * pow(2 * PI, 6));

    // q -> k1 k2 phase space
    dbl_type e1 = (q2 + m1 * m1 - m2 * m2) / (2 * sqrt(q2));
    dbl_type e2 = (q2 + m2 * m2 - m1 * m1) / (2 * sqrt(q2));
    mom = sqrt(e1 * e1 - m1 * m1);
    cosT = random->rand(-1, 1);
    sinT = sqrt(1. - cosT * cosT);
    WT *= 2;
    phi = random->rand(0, 2 * PI);
    sinPhi = sin(phi);
    cosPhi = cos(phi);
    WT *= 2 * PI;

    k1[3] = e1;
    k1[0] = mom * sinT * cosPhi;
    k1[1] = mom * sinT * sinPhi;
    k1[2] = mom * cosT;

    k2[3] = e2;
    k2[0] = -mom * sinT * cosPhi;
    k2[1] = -mom * sinT * sinPhi;
    k2[2] = -mom * cosT;

    WT *= mom / (sqrt(q2) * 4 * pow(2 * PI, 6));

    // boost
    apply_boost_to(q, k1);
    apply_boost_to(q, k2);

    WT *= pow(2 * PI, 3+4);
    return WT;
}
