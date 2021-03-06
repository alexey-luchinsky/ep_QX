/* 
 * File:   RamboEP.c++
 * Author: luchinsky
 * 
 * Created on January 26, 2017, 12:21 PM
 */

#include "RamboEP.h"

RamboEP::RamboEP(dbl_type ecm_, Random *random_, dbl_type _m1, dbl_type _m2, dbl_type _m3) {
    ecm=ecm_;
    m1=_m1;
    m2=_m2;
    m3=_m3;
    random=random_;
    ram2 = new Rambo2(m1,m2,random);
    ram3 = new Rambo3(m1,m2,m3,random);
    S=ecm*ecm;
    set_v4(P,0,0,-ecm/2,ecm/2);
    set_v4(kIn,0,0,ecm/2,ecm/2);
    minQ2=0; maxQ2=S;
    minY=0; maxY=1;
    PI=acos(-1);
    alpha = 1./137;
    nFault=0;
}

RamboEP::~RamboEP() {
    delete ram2;
    delete ram3;
}

void RamboEP::set_minQ2(dbl_type _minQ2) {
    if(_minQ2<0) {
        cout<<"RamboEP::set_minQ2: _minQ2="<<_minQ2<<"<0, skipping"<<endl;
        return;
    };
    minQ2=_minQ2;
}

void RamboEP::set_maxQ2(dbl_type _maxQ2) {
    if(_maxQ2>S) {
        cout<<"RamboEP::set_maxQ2: _minQ2="<<_maxQ2<<">S="<<S<<", skipping"<<endl;
        return;
    };
    maxQ2=_maxQ2;
}

void RamboEP::set_minmaxQ2(dbl_type _minQ2, dbl_type _maxQ2) {
    set_minQ2(_minQ2); set_maxQ2(_maxQ2);
}

bool RamboEP::next(dbl_type x) {
//    cout<<"RamboEP::next, start"<<endl;
    set_v4(Pg,P);  mult(Pg,x);
    wt=1./(16*PI*PI);
    if(maxY<minY) {
        wt=0; wT=0; wL=0; return false;        
    };
    Y = random->rand(minY, maxY);
    wt *= maxY-minY;
    dbl_type _maxQ2=min(maxQ2,S*x*Y);
    if(_maxQ2<minQ2) {
        wt=0; wT=0; wL=0; return false;
    };
    Q2 = random->rand(minQ2, _maxQ2);
    wt *= _maxQ2-minQ2;
    W2 = -Q2+x*S*Y;
    dbl_type omega = ecm*(1-Y+Q2/S)/2;
    dbl_type cosT = (1-Y-Q2/S)/(1-Y+Q2/S), sinT=sqrt(1-cosT*cosT);
    set_v4(kOut, 0, omega*sinT, omega*cosT, omega);
    subtract(kIn, kOut, q);
    wT = wt*8*PI*alpha/Q2*(1+pow(1-Y,2))/pow(Y,2);
    wL = -wt*32*PI*alpha/Q2*(1-Y)/pow(Y,2);
    return true;
}

bool RamboEP::next(dbl_type (&kp_)[4], dbl_type (&k1)[4], dbl_type (&k2)[4], dbl_type x) {
    dbl_type Phard[4], _k1[4], _k2[4];
    bool ok=true;
    if(!next(x)) ok=false;    
    if(sqrt(W2)<m1+m2) ok=false;

    set_v4(kp_, kOut);
    sum(q,Pg,Phard);      // Phard = q+p
    dbl_type _wt=ram2->next(sqrt(W2),_k1, _k2);
    if(!(_wt>=0)) ok=false;
    if(!ok) {
        nFault++;
        wt=0; wT=0; wL=0;
        return false;
    };
    set_v4(k1,_k1); set_v4(k2,_k2);
    apply_boost_to(Phard,k1);  apply_boost_to(Phard, k2);
    if(!are_equal( mass2(k1), m1*m1)) {
        cout<<" Error at RamEP::next\t k1^2 != m1^2"<<endl;
        cout<<" W="<<sqrt(W2)<<endl;
        println_4v("Phard",Phard);
        println_4v("_k1",_k1);
        println_4v("k1",k1);
    };
    assert(are_equal( mass2(k1), m1*m1));
    assert(are_equal( mass2(k2), m2*m2));
    assert(are_equal( sum_mass2(k1,k2),W2));
    wt *= _wt;
    wT *= _wt;
    wL *= _wt;
    return true;
}



