#include "kinematics/RamboEP.h"
#include "kinematics/algebra.h"
using namespace std;

dbl_type M = 3.1, M2 = M*M;

inline void assert_is_zero(dbl_type a) { assert(is_zero(a));}
inline void assert_are_equal(dbl_type a, dbl_type b) { assert(are_equal(a,b));}

int main(void) {
    cout << "3S1_1c" << endl;
    Random *random_generator = new Random();
    dbl_type ecm = 50, S = ecm*ecm;
    dbl_type x = 0.5;
    dbl_type P[4], PgOut[4], k[4], kp[4], q[4], pPsi[4], kg[4];
    RamboEP ramEP(ecm, random_generator, M, 0);
    set_v4(P, ramEP.P);
    set_v4(PgOut, P);
    mult(PgOut, 1 - x);
    set_v4(k, ramEP.kIn);
    set_v4(kp, ramEP.kOut);
    assert_is_zero(mass2(P));
    assert_is_zero(mass2(k));

    ramEP.next(kp, pPsi, kg, x);
    // check kinematics
    assert_is_zero(mass2(kp));
    assert_is_zero(mass2(kg));
    assert_are_equal(mass2(pPsi),M2);
    dbl_type PtotIn[4], PtotOut[4];
    sum(P, k, PtotIn);
    sum(kp, pPsi, kg, PtotOut); sum(PtotOut, PgOut, PtotOut);
    assert_is_zero(PtotIn[0]);  assert_is_zero(PtotOut[0]);
    assert_is_zero(PtotIn[1]);  assert_is_zero(PtotOut[1]);
    assert_is_zero(PtotIn[2]);  assert_is_zero(PtotOut[2]);
    assert_are_equal(PtotIn[3], ecm); assert_are_equal(PtotOut[3], ecm);




    return 0;
}
