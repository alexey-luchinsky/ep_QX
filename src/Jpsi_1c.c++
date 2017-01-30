#include "kinematics/RamboEP.h"
#include "kinematics/algebra.h"
using namespace std;

dbl_type M = 3.1, M2 = M*M;
dbl_type PI = acos(-1), alpha = 1. / 137, alphas = 0.3;

inline void assert_is_zero(dbl_type a) {
    assert(is_zero(a));
}

inline void assert_are_equal(dbl_type a, dbl_type b) {
    assert(are_equal(a, b));
}

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

    for (int iEv = 0; iEv < 10; ++iEv) {
        ramEP.next(kp, pPsi, kg, x);
        subtract(ramEP.kIn, ramEP.kOut, q);
        // check kinematics
        assert_is_zero(mass2(kp));  assert_is_zero(mass2(kg)); assert_are_equal(mass2(pPsi), M2);
        dbl_type PtotIn[4], PtotOut[4];
        sum(P, k, PtotIn); sum(kp, pPsi, kg, PtotOut);  sum(PtotOut, PgOut, PtotOut);
        assert_is_zero(PtotIn[0]);
        assert_is_zero(PtotOut[0]);
        assert_is_zero(PtotIn[1]);
        assert_is_zero(PtotOut[1]);
        assert_is_zero(PtotIn[2]);
        assert_is_zero(PtotOut[2]);
        assert_are_equal(PtotIn[3], ecm);
        assert_are_equal(PtotOut[3], ecm);

        // gamma g -> JP g Mandelstam variables
        dbl_type Pg[4];
        set_v4(Pg, ramEP.Pg);
        dbl_type Q2 = -mass2(q), s = sum_mass2(Pg, q), t = subtract_mass2(q, kg), u = subtract_mass2(q, pPsi);
        assert_are_equal(s + t + u, M2 - Q2);


        dbl_type matr2T = (-2097152 * alpha * pow(alphas, 2) * pow(M, 2) * pow(PI, 3) * pow(s, -2) * pow(Q2 + s, -2) * pow(u, -2)*
                ((s + u) * pow(Q2, 6) + (s + u) * pow(s, 4) * pow(u, 2) + u * pow(Q2, 4)*(-5 * s * u - 8 * pow(s, 2) + pow(u, 2)) + Q2 * u * pow(s, 3)*(4 * s * u + pow(s, 2) + 3 * pow(u, 2)) -
                pow(Q2, 5)*(5 * s * u + pow(s, 2) + 3 * pow(u, 2)) + s * pow(Q2, 3)*(-4 * u * pow(s, 2) + pow(s, 3) - s * pow(u, 2) + 3 * pow(u, 3)) +
                pow(Q2, 2) * pow(s, 2)*(u * pow(s, 2) + pow(s, 3) + 4 * s * pow(u, 2) + 4 * pow(u, 3))) * pow(s + u, -1) * pow(2 * Q2 + s + u, -2)) / 9.;

        dbl_type matr2L = (-1048576 * alpha * Q2 * pow(alphas, 2) * pow(M, 2) * pow(PI, 3) * pow(s, -2) * pow(Q2 + s, -2) * pow(u, -2) * pow(s + u, -1)*
                (-2 * Q2 * u * (s + 2 * u) * pow(s, 3) - pow(s, 4) * pow(u, 2) - 2 * s * pow(Q2, 3)*(3 * s * u + 2 * pow(s, 2) + pow(u, 2)) - pow(Q2, 2) * pow(s, 2)*(8 * s * u + 2 * pow(s, 2) + 7 * pow(u, 2)) +
                2 * pow(Q2, 4) * pow(s + u, 2)) * pow(2 * Q2 + s + u, -2)) / 9.;

        dbl_type matr2 = matr2L + ramEP.wL + matr2T * ramEP.wT;
        cout << "matr2=" << matr2 << endl;
        assert(matr2 > 0);
    };
    return 0;
}
