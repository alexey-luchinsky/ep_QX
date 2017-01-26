/* 
 * File:   rambo3.h
 * Author: luchinsky
 *
 * Created on March 21, 2016, 11:34 AM
 * 3D LIPS q k3 -> k1 k2 k3
 */

#ifndef RAMBO3_H
#define RAMBO3_H

#include <iostream>
#include <stdlib.h>
#include "algebra.h"
#include "Random.h"

using namespace std;

class Rambo3 {
private:
    Random *random;
    dbl_type m1, m2, m3;
    dbl_type q[4];
    dbl_type PI;
    dbl_type min_cos, max_cos;

public:
    Rambo3() { }

    Rambo3(dbl_type m1_, dbl_type m2_, dbl_type m3_, Random *random_);

    dbl_type next(dbl_type ecm, dbl_type (&k1)[4], dbl_type (&k2)[4], dbl_type (&k3)[4]);
};

#endif	/* RAMBO3_H */

