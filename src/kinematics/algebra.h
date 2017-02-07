#ifndef GGGPSIPSI_ALGEBRA_H
#define GGGPSIPSI_ALGEBRA_H

#include <iostream>
#include <complex>
#include <math.h>
#include <stdexcept>
//#include "color_util.h"
#include <sstream>

//#define NDEBUG
#include <assert.h>
//#include "common_util.h"

using namespace std;

/************************** DATA TYPES **************************/

//#define USE_EXT_DOUBLE

#ifdef USE_EXT_DOUBLE
typedef long double dbl_type;
#else
typedef double dbl_type;
#endif


typedef complex<dbl_type> cmplx;

#define IMAGE_ONE cmplx((dbl_type) 0.0, (dbl_type) 1.0);
#define COMPLEX_ZERO cmplx((dbl_type) 0.0, (dbl_type) 0.0);
#define COMPLEX_ONE cmplx((dbl_type) 1.0, (dbl_type) 0.0);


inline
cmplx operator+(const cmplx &v1, const int &v2) {
    return v1 + ((dbl_type) v2);
}

inline
cmplx operator+(const int &v2, const cmplx &v1) {
    return v1 + ((dbl_type) v2);
}

inline
cmplx operator*(const cmplx &v1, const int &v2) {
    return v1 * ((dbl_type) v2);
}

inline
cmplx operator*(const int &v2, const cmplx &v1) {
    return v1 * ((dbl_type) v2);
}

inline
cmplx operator/(const cmplx &v1, const int &v2) {
    return v1 / ((dbl_type) v2);
}

inline
cmplx operator/(const int &v2, const cmplx &v1) {
    cmplx tmp = COMPLEX_ONE;
    tmp /= v1;
    tmp *= (dbl_type) v2;
    return tmp;
}

/************************** VECTOR ALGEBRA **************************/

//don't inline these functions:
dbl_type sp(dbl_type (&v1)[4], dbl_type (&v2)[4]);

cmplx sp(cmplx (&v1)[4], cmplx (&v2)[4]);

cmplx sp(cmplx (&v1)[4], dbl_type (&v2)[4]);

cmplx sp(dbl_type (&v2)[4], cmplx (&v1)[4]);

template<typename T>
inline
T mass2(T (&v4)[4]) {
    return sp(v4, v4);
}

template<typename T>
inline
T mass(T (&v4)[4]) {
    return sqrt(mass2(v4));
}

inline
void set_v4(dbl_type (&v)[4], dbl_type x, dbl_type y, dbl_type z, dbl_type t) {
    v[0] = x;
    v[1] = y;
    v[2] = z;
    v[3] = t;
}

inline 
void set_v4(dbl_type (&vOut)[4], dbl_type vIn[4]) {
    vOut[0]=vIn[0];
    vOut[1]=vIn[1];
    vOut[2]=vIn[2];
    vOut[3]=vIn[3];
}

void subtract(dbl_type (&v1)[4], dbl_type (&v2)[4], dbl_type (&answ)[4]);

void subtract(cmplx (&v1)[4], cmplx (&v2)[4], cmplx (&answ)[4]);

void sum(dbl_type (&v1)[4], dbl_type (&v2)[4], dbl_type (&answ)[4]);

void sum(dbl_type (&v1)[4], dbl_type (&v2)[4], dbl_type (&v3)[4], dbl_type (&answ)[4]);

void sum(cmplx (&v1)[4], cmplx (&v2)[4], cmplx (&answ)[4]);

void mult(dbl_type (&v1)[4], dbl_type factor);

dbl_type subtract_mass2(dbl_type (&v1)[4], dbl_type (&v2)[4]);

dbl_type sum_mass2(dbl_type (&v1)[4], dbl_type (&v2)[4]);

dbl_type sum_mass2(dbl_type (&v1)[4], dbl_type (&v2)[4], dbl_type (&v3)[4]);

dbl_type angle(dbl_type (&v1)[4], dbl_type (&v2)[4]);

inline
dbl_type sum_mass(dbl_type (&v1)[4], dbl_type (&v2)[4]) { return sqrt(sum_mass2(v1, v2)); }

dbl_type pT_squared(dbl_type (&v)[4]);

inline
dbl_type pT(dbl_type (&v)[4]) { return sqrt(pT_squared(v)); }

inline
dbl_type rapidity(dbl_type (&v)[4]) { return log((v[3] - v[2]) / (v[3] + v[2])) / 2.0; }

inline
dbl_type sum_rapidity(dbl_type (&v1)[4], dbl_type (&v2)[4]) {
    return log((v1[3] + v2[3] - v1[2] - v2[2]) / (v1[3] + v2[3] + v1[2] + v2[2])) / 2.0;
}

dbl_type sum_pT_squared(dbl_type (&v1)[4], dbl_type (&v2)[4]);

inline
dbl_type sum_pT(dbl_type (&v1)[4], dbl_type (&v2)[4]) { return sqrt(sum_pT_squared(v1, v2)); }

inline
dbl_type __atan2___(dbl_type y, dbl_type x) {
    if (x != 0) return atan2(y, x);
    if (y == 0) return 0;
    if (y > 0) return M_PI / 2;
    else return -M_PI / 2;
}

inline
dbl_type phi(dbl_type (&v)[4]) {
    return v[0] == 0.0 && v[1] == 0.0 ? 0.0 : __atan2___(v[1], v[0]);
}

/************************** BOOST **************************/

void apply_boost_to(dbl_type bx, dbl_type by, dbl_type bz, cmplx (&answ)[4]);

void apply_boost_to(dbl_type bx, dbl_type by, dbl_type bz, dbl_type (&answ)[4]);

void apply_boost_to(dbl_type (&boost)[4], cmplx (&answ)[4]);

void apply_boost_to(dbl_type (&boost)[4], dbl_type (&answ)[4]);


/************************** POLARIZATIONS **************************/

void set_polarization_to_momentum(cmplx (&pol)[4], dbl_type (&mom)[4]);

void rotate_euler(dbl_type phi, dbl_type theta, dbl_type ksi, cmplx (&answ)[4]);

void set_gluon_polarization(int state, dbl_type (&momentum)[4], cmplx (&answ)[4]);

void set_gluon_polarization(int state, dbl_type (&momentum)[4], cmplx (&answ)[4], cmplx (cans)[4]);

void set_psi_polarization(int state, dbl_type (&momentum)[4], cmplx (&answ)[4]);

void set_psi_polarization(int state, dbl_type (&momentum)[4], cmplx (&answ)[4], cmplx (&cans)[4]);


/************************** UTIL **************************/

template<typename T>
inline
void println_4v(const string &name, T (&v4)[4]) {
    cout << name << ": ";
    cout << "{" << v4[3] << ", " << v4[0] << ", " << v4[1] << ", " << v4[2] << "}";
    cout << ",  M2 = " << mass2(v4) << ",  " << " M = " << mass(v4) << endl;
}

template<typename T>
inline
void println_4v_math(const string &name, T (&v4)[4]) {
    cout << name << "= ";
    cout << "{" << v4[3] << ", " << v4[0] << ", " << v4[1] << ", " << v4[2] << "};"<<endl;
}


template<typename T>
inline
void println_4v(T (&v4)[4]) {
    println_4v("", v4);
}


/************************** ASSERT **************************/

//assertions
#define DEFAULT_ZERO_PRECISION 1e-3

inline
bool is_zero(dbl_type a, double precision) {
    return a < ((dbl_type) precision) && a > -((dbl_type) precision);
}

inline
bool is_zero(cmplx a, double precision) {
    return is_zero(a.real(), precision) && is_zero(a.imag(), precision);
}

inline
bool is_zero(dbl_type a) {
    return is_zero(a, DEFAULT_ZERO_PRECISION);
}

inline
bool is_zero(cmplx a) {
    return is_zero(a, DEFAULT_ZERO_PRECISION);
}

inline
bool is_real(cmplx a) {
    return is_zero(a.imag());
}

inline
bool are_equal(dbl_type a, dbl_type b, double precision) {
    return is_zero(a - b, precision);
}

inline
bool are_equal(dbl_type a, dbl_type b) {
    return are_equal(a, b, DEFAULT_ZERO_PRECISION);
}

#endif