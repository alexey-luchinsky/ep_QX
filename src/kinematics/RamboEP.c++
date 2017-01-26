/* 
 * File:   RamboEP.c++
 * Author: luchinsky
 * 
 * Created on January 26, 2017, 12:21 PM
 */

#include "RamboEP.h"

RamboEP::RamboEP(dbl_type ecm_, Random *random_) {
    ecm=ecm_;
    S=ecm*ecm;
    random=random_;
    set_v4(P,0,0,-ecm/2,ecm/2);
    set_v4(kIn,0,0,ecm/2,ecm/2);
    minQ2=0; maxQ2=S;
    minY=0; maxY=1;
    PI=acos(-1);
}

void RamboEP::next() {
    dbl_type WT=1;
    Y = random->rand(minY, maxY); WT *= maxY-minY;
    Q2 = random->rand(minQ2, maxQ2); WT *= maxQ2 - minQ2;
    WT *= 1./(16*PI*PI);
    dbl_type omega = (Q2+S*(1-Y))/(2*ecm);
    dbl_type cosT = (S*(1-Y)-Q2)/(S*(1-Y)+Q2), sinT=sqrt(1-cosT*cosT);
    set_v4(kOut, 0, omega*sinT, omega*cosT, omega);
    subtract(kIn, kOut, q);
    wT = (1+pow(1-Y,2))/Y; wL=-4*(1-Y)/Y;
}



