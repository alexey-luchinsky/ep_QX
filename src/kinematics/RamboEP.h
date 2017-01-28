/* 
 * File:   RamboEP.h
 * Author: luchinsky
 *
 * Created on January 26, 2017, 12:21 PM
 */

#ifndef RAMBOEP_H
#define	RAMBOEP_H

#include "Random.h"
#include "Rambo2.h"
#include "Rambo3.h"


class RamboEP {
public:
    RamboEP(dbl_type _ecm, Random *random_, dbl_type _m1=0, dbl_type _m2=0, dbl_type _m3=0);
    ~RamboEP();
    bool next(dbl_type x=1);
    bool next(dbl_type (&_kp)[4], dbl_type (&k1)[4], dbl_type (&k2)[4], dbl_type x=1);
//private:
    Random *random;
    Rambo2 *ram2;
    Rambo3 *ram3;
    dbl_type ecm, S;
    dbl_type m1, m2, m3;
    dbl_type Q2, minQ2, maxQ2;
    dbl_type Y, minY, maxY;
    dbl_type P[4], kIn[4], kOut[4], q[4], Pg[4], PhardG[4];
    dbl_type W2;
    dbl_type wT, wL, wt;
    dbl_type PI;
    dbl_type alpha;
};

#endif	/* RAMBOEP_H */

