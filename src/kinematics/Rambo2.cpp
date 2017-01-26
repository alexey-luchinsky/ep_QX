/* 
 * File:   Rambo2.cpp
 * Author: luchinsky
 * 
 * Created on January 26, 2017, 11:00 AM
 */

#include "Rambo2.h"

Rambo2::Rambo2(dbl_type m1_, dbl_type m2_, Random *random_) {
    m1 = m1_;
    m2 = m2_;
    random = random_;
    PI = acos(-1);

}

dbl_type Rambo2::next(dbl_type ecm, dbl_type(&k1)[4], dbl_type(&k2)[4]) {
    dbl_type WT = 1;
    dbl_type e1 = (ecm * ecm + m1 * m1 - m2 * m2) / (2 * ecm), e2 = (ecm * ecm + m2 * m2 - m1 * m1) / (2 * ecm);
    dbl_type mom = sqrt(e1*e1-m1*m1);
    k1[3] = e1; k2[3]=e2;
    dbl_type cosTH = random->rand(-1, 1), sinTH = sqrt(1 - cosTH * cosTH);
    WT *= 2;
    dbl_type phi=random->rand(0,2*PI);
    WT *= 2*PI;
    phi=random->rand(0,2*PI);
    k1[0]=mom*sinTH*cos(phi);    k2[0] = -k1[0];
    k1[1]=mom*sinTH*sin(phi);    k2[1]=-k1[1];
    k1[2]=mom*cosTH;    k2[2]=-k1[2];
    WT *= mom/ecm/(16*PI*PI);
    return WT;
}


