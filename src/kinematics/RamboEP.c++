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
    alpha = 1./137;
}

void RamboEP::next() {
    wt=1./(16*PI*PI);
    Y = random->rand(minY, maxY);
    wt *= maxY-minY;
    dbl_type _maxQ2=min(maxQ2,S*Y);
    Q2 = random->rand(minQ2, _maxQ2);
    wt *= _maxQ2-minQ2;
    W2 = -Q2+S*Y;
    dbl_type omega = (Q2+S*(1-Y))/(2*ecm);
    dbl_type cosT = (S*(1-Y)-Q2)/(S*(1-Y)+Q2), sinT=sqrt(1-cosT*cosT);
    set_v4(kOut, 0, omega*sinT, omega*cosT, omega);
    subtract(kIn, kOut, q);
    wT = alpha/(pow(2*PI,2)*Q2)*(1+pow(1-Y,2))/pow(Y,2)*(maxY-minY)*(_maxQ2-minQ2);
    wL = -2*alpha/(2*PI*PI*Q2)*(1-Y)/pow(Y,2)*(maxY-minY)*(_maxQ2-minQ2);
}



