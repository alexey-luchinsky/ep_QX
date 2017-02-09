#include "algebra.h"


dbl_type sp(dbl_type (&v1)[4], dbl_type (&v2)[4]) {
    return v1[3] * v2[3] - v1[0] * v2[0] - v1[1] * v2[1] - v1[2] * v2[2];
}

cmplx sp(cmplx (&v1)[4], cmplx (&v2)[4]) {
    return v1[3] * v2[3] - v1[0] * v2[0] - v1[1] * v2[1] - v1[2] * v2[2];
}

cmplx sp(cmplx (&v1)[4], dbl_type (&v2)[4]) {
    return v1[3] * v2[3] - v1[0] * v2[0] - v1[1] * v2[1] - v1[2] * v2[2];
}

cmplx sp(dbl_type (&v2)[4], cmplx (&v1)[4]) {
    return sp(v1, v2);
}

void subtract(dbl_type (&v1)[4], dbl_type (&v2)[4], dbl_type (&answ)[4]) {
    for (int i = 0; i < 4; ++i)
        answ[i] = v1[i] - v2[i];
}

void subtract(cmplx (&v1)[4], cmplx (&v2)[4], cmplx (&answ)[4]) {
    for (int i = 0; i < 4; ++i)
        answ[i] = v1[i] - v2[i];
}

void sum(dbl_type (&v1)[4], dbl_type (&v2)[4], dbl_type (&answ)[4]) {
    for (int i = 0; i < 4; ++i)
        answ[i] = v1[i] + v2[i];
}

void sum(dbl_type (&v1)[4], dbl_type (&v2)[4], dbl_type (&v3)[4], dbl_type (&answ)[4]) {
    for (int i = 0; i < 4; ++i)
        answ[i] = v1[i] + v2[i] + v3[i];
}

void sum(cmplx (&v1)[4], cmplx (&v2)[4], cmplx (&answ)[4]) {
    for (int i = 0; i < 4; ++i)
        answ[i] = v1[i] + v2[i];
}

void mult(dbl_type (&v1)[4], dbl_type factor) {
    for (int i = 0; i < 4; ++i)
        v1[i] *= factor;
}

dbl_type subtract_mass2(dbl_type (&v1)[4], dbl_type (&v2)[4]) {
    dbl_type e = v1[3] - v2[3];
    dbl_type x = v1[0] - v2[0];
    dbl_type y = v1[1] - v2[1];
    dbl_type z = v1[2] - v2[2];
    return e * e - x * x - y * y - z * z;
}

dbl_type sum_mass2(dbl_type (&v1)[4], dbl_type (&v2)[4]) {
    dbl_type e = v1[3] + v2[3];
    dbl_type x = v1[0] + v2[0];
    dbl_type y = v1[1] + v2[1];
    dbl_type z = v1[2] + v2[2];
    return e * e - x * x - y * y - z * z;
}

dbl_type sum_mass2(dbl_type (&v1)[4], dbl_type (&v2)[4], dbl_type (&v3)[4]) {
    dbl_type e = v1[3] + v2[3] + v3[3];
    dbl_type x = v1[0] + v2[0] + v3[0];
    dbl_type y = v1[1] + v2[1] + v3[1];
    dbl_type z = v1[2] + v2[2] + v3[2];
    return e * e - x * x - y * y - z * z;
}

dbl_type pT_squared(dbl_type (&v)[4]) {
    return v[0] * v[0] + v[1] * v[1];
}

dbl_type sum_pT_squared(dbl_type (&v1)[4], dbl_type (&v2)[4]) {
    dbl_type px = v1[0] + v2[0];
    dbl_type py = v1[1] + v2[1];
    return px * px + py * py;
}

dbl_type angle(dbl_type (&v1)[4], dbl_type (&v2)[4]) {
    dbl_type sp = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    if (sp <= 0)
        return 0.0;
    dbl_type m1 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
    dbl_type m2 = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];
    dbl_type cos_ = sp / sqrt(m1 * m2);
    if (cos_ > 1.0)
        cos_ = 1.0;
    if (cos_ < -1.0)
        cos_ = -1.0;
    return acos(cos_);
}

void set_polarization_to_momentum(cmplx (&pol)[4], dbl_type (&mom)[4]) {
    for (int i = 0; i < 4; ++i)
        pol[i] = cmplx(mom[i], (dbl_type) 0.0);
}


void rotate_euler(dbl_type phi, dbl_type theta, dbl_type ksi, cmplx (&answ)[4]) {
    dbl_type sp = sin(phi);
    dbl_type st = sin(theta);
    dbl_type sk = sin(ksi);
    dbl_type cp = cos(phi);
    dbl_type ct = cos(theta);
    dbl_type ck = cos(ksi);

    cmplx x = (ck * ct * cp - sk * sp) * answ[0] + (-sk * ct * cp - ck * sp) * answ[1] + st * cp * answ[2];
    cmplx y = (ck * ct * sp + sk * cp) * answ[0] + (-sk * ct * sp + ck * cp) * answ[1] + st * sp * answ[2];
    cmplx z = -ck * st * answ[0] + sk * st * answ[1] + ct * answ[2];

    answ[0] = x;
    answ[1] = y;
    answ[2] = z;
}

void set_gluon_polarization(int state, dbl_type (&momentum)[4], cmplx (&answ)[4]) {
    cmplx tmp_eps1[4];
    cmplx tmp_eps2[4];
    dbl_type sqrt2 = sqrt((dbl_type) 2.0);
    dbl_type zero = (dbl_type) 0;

    tmp_eps1[3] = COMPLEX_ZERO;
    tmp_eps1[0] = cmplx(-1.0 / sqrt2, zero);
    tmp_eps1[1] = cmplx(zero, -1.0 / sqrt2);
    tmp_eps1[2] = COMPLEX_ZERO;

    tmp_eps2[3] = COMPLEX_ZERO;
    tmp_eps2[0] = cmplx(1.0 / sqrt2, zero);
    tmp_eps2[1] = cmplx(zero, -1.0 / sqrt2);
    tmp_eps2[2] = COMPLEX_ZERO;

    dbl_type px = momentum[0];
    dbl_type py = momentum[1];
    dbl_type pz = momentum[2];
    dbl_type phi = atan2(py, px);
    dbl_type theta = acos(pz / sqrt(px * px + py * py + pz * pz));

    rotate_euler(phi, theta, -phi, tmp_eps1);
    rotate_euler(phi, theta, -phi, tmp_eps2);

    switch (state) {
        case 0:
            copy(tmp_eps1, tmp_eps1 + 4, answ);
            break;
        case 1:
            copy(tmp_eps2, tmp_eps2 + 4, answ);
            break;
        default:
            throw runtime_error("unknown polarization");
    }
}

void set_gluon_polarization(int state, dbl_type (&momentum)[4], cmplx (&answ)[4], cmplx (cans)[4]) {
    set_gluon_polarization(state, momentum, answ);
    for(int i=0; i<4; ++i) cans[i]=conj(answ[i]);
}

void lepton_current(int state, dbl_type (&p1)[4], dbl_type (&p2)[4], cmplx (&eps)[4]) {
    const cmplx II=IMAGE_ONE;
    cmplx norm=-4*sp(p1,p2);;
    if(state==0) {
        eps[0]=-((-p1[2] + p1[3])*(-p2[0] - II*p2[1])) - (-p1[0] + II*p1[1])*(-p2[2] + p2[3]);
        eps[1]=II*(-p1[2] + p1[3])*(-p2[0] - II*p2[1]) - II*(-p1[0] + II*p1[1])*(-p2[2] + p2[3]);
        eps[2]=-((p1[0] - II*p1[1])*(-p2[0] - II*p2[1])) - (-p1[2] + p1[3])*(-p2[2] + p2[3]);
        eps[3]=(-p1[0] + II*p1[1])*(-p2[0] - II*p2[1]) + (-p1[2] + p1[3])*(-p2[2] + p2[3]);
        norm=norm/(4*(p1[2] - p1[3])*(p2[2] - p2[3])*(p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2] - p1[3]*p2[3]));
    }
    else if(state==1) {
	 eps[0]=(p1[2] + p1[3])*(p2[0] + II*p2[1]) + (p1[0] - II*p1[1])*(p2[2] + p2[3]);
	 eps[1]=-(II*(p1[2] + p1[3])*(p2[0] + II*p2[1])) + II*(p1[0] - II*p1[1])*(p2[2] + p2[3]);
	 eps[2]=(-p1[0] + II*p1[1])*(p2[0] + II*p2[1]) + (p1[2] + p1[3])*(p2[2] + p2[3]);
	 eps[3]=(p1[0] - II*p1[1])*(p2[0] + II*p2[1]) + (p1[2] + p1[3])*(p2[2] + p2[3]);
         
        norm=norm/(4*(p1[2] + p1[3])*(p2[2] + p2[3])*(p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2] - p1[3]*p2[3]));
    }
    else {
        cout<<"Wrong lepton state "<<state<<endl;
        return;
    };
    for(int i=0; i<4; ++i) eps[i] = sqrt(norm)*eps[i];
}


void apply_boost_to(dbl_type bx, dbl_type by, dbl_type bz, cmplx (&answ)[4]) {
    dbl_type bxx = bx * bx;
    dbl_type byy = by * by;
    dbl_type bzz = bz * bz;

    dbl_type b2 = bxx + byy + bzz;
    if (b2 == (dbl_type) 0.0)
        return;
    assert(b2 <= 1.0);

    dbl_type gamma = 1.0 / sqrt(1 - b2);
    dbl_type gb2 = (gamma - 1.0) / b2;
    dbl_type gb2xy = gb2 * bx * by;
    dbl_type gb2xz = gb2 * bx * bz;
    dbl_type gb2yz = gb2 * by * bz;

    dbl_type gbx = gamma * bx;
    dbl_type gby = gamma * by;
    dbl_type gbz = gamma * bz;

    cmplx e2 = answ[3];
    cmplx px2 = answ[0];
    cmplx py2 = answ[1];
    cmplx pz2 = answ[2];

    answ[3] = gamma * e2 + gbx * px2 + gby * py2 + gbz * pz2;
    answ[0] = gbx * e2 + gb2 * bxx * px2 + px2 + gb2xy * py2 + gb2xz * pz2;
    answ[1] = gby * e2 + gb2 * byy * py2 + py2 + gb2xy * px2 + gb2yz * pz2;
    answ[2] = gbz * e2 + gb2 * bzz * pz2 + pz2 + gb2yz * py2 + gb2xz * px2;
}

void apply_boost_to(dbl_type (&boost)[4], cmplx (&answ)[4]) {
    dbl_type bx = (boost[0] / boost[3]);
    dbl_type by = (boost[1] / boost[3]);
    dbl_type bz = (boost[2] / boost[3]);

    apply_boost_to(bx, by, bz, answ);
}

void apply_boost_to(dbl_type bx, dbl_type by, dbl_type bz, dbl_type (&answ)[4]) {
    dbl_type bxx = bx * bx;
    dbl_type byy = by * by;
    dbl_type bzz = bz * bz;

    dbl_type b2 = bxx + byy + bzz;
    if (b2 == (dbl_type) 0.0)
        return;
    assert(b2 <= 1.0);

    dbl_type gamma = 1.0 / sqrt(1 - b2);
    dbl_type gb2 = (gamma - 1.0) / b2;
    dbl_type gb2xy = gb2 * bx * by;
    dbl_type gb2xz = gb2 * bx * bz;
    dbl_type gb2yz = gb2 * by * bz;

    dbl_type gbx = gamma * bx;
    dbl_type gby = gamma * by;
    dbl_type gbz = gamma * bz;

    dbl_type e2 = answ[3];
    dbl_type px2 = answ[0];
    dbl_type py2 = answ[1];
    dbl_type pz2 = answ[2];

    answ[3] = gamma * e2 + gbx * px2 + gby * py2 + gbz * pz2;
    answ[0] = gbx * e2 + gb2 * bxx * px2 + px2 + gb2xy * py2 + gb2xz * pz2;
    answ[1] = gby * e2 + gb2 * byy * py2 + py2 + gb2xy * px2 + gb2yz * pz2;
    answ[2] = gbz * e2 + gb2 * bzz * pz2 + pz2 + gb2yz * py2 + gb2xz * px2;
}

void apply_boost_to(dbl_type (&boost)[4], dbl_type (&answ)[4]) {
    dbl_type bx = (dbl_type) (boost[0] / boost[3]);
    dbl_type by = (dbl_type) (boost[1] / boost[3]);
    dbl_type bz = (dbl_type) (boost[2] / boost[3]);
    apply_boost_to(bx, by, bz, answ);
}

void set_psi_polarization(int state, dbl_type (&momentum)[4], cmplx (&answ)[4]) {
    for (int i = 0; i < 4; ++i)
        answ[i] = COMPLEX_ZERO;
    answ[state] = COMPLEX_ONE;
    apply_boost_to(momentum, answ);
}

void set_psi_polarization(int state, dbl_type (&momentum)[4], cmplx (&answ)[4], cmplx (&cans)[4]) {
    set_psi_polarization(state, momentum, answ);
    for(int i=0; i<4; ++i) cans[i]=conj(answ[i]);
}