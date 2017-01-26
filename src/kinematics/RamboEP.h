/* 
 * File:   RamboEP.h
 * Author: luchinsky
 *
 * Created on January 26, 2017, 12:21 PM
 */

#ifndef RAMBOEP_H
#define	RAMBOEP_H

#include "Random.h"


class RamboEP {
public:
    RamboEP(dbl_type _ecm, Random *random_);
    void next();
//private:
    Random *random;
    dbl_type ecm, S;
    dbl_type Q2, minQ2, maxQ2;
    dbl_type Y, minY, maxY;
    dbl_type P[4], kIn[4], kOut[4], q[4];
    dbl_type W2;
    dbl_type wT, wL;
    dbl_type PI;
    dbl_type alpha = 1./137;
};

#endif	/* RAMBOEP_H */

