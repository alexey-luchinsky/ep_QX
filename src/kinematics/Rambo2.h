/* 
 * File:   Rambo2.h
 * Author: luchinsky
 *
 * Created on January 26, 2017, 11:00 AM
 */

#ifndef RAMBO2_H
#define	RAMBO2_H

#include "Random.h"
#include "algebra.h"

class Rambo2 {
public:
    Rambo2() { };
    Rambo2(dbl_type m1_, dbl_type m2_, Random *random_);
    dbl_type next(dbl_type ecm, dbl_type (&k1)[4], dbl_type (&k2)[4]);
private:
    Random *random;
    dbl_type m1, m2;
    dbl_type PI;
};

#endif	/* RAMBO2_H */

