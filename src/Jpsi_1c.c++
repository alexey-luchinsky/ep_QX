#include "kinematics/RamboEP.h"

int main(void) {
    dbl_type ecm = 10, M = 3.1, mc = M / 2, M2 = M*M;
    dbl_type x = 1;
    Random *random_generator = new Random();
    RamboEP ramEP(ecm, random_generator, M, 0);
    dbl_type k[4], P[4], kp[4], k1[4], k2[4], p[4], k3[4];
    set_v4(k, ramEP.kIn);
    set_v4(P, ramEP.P);

    int nEv=100;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        cout << "(*====== Event " << iEv << " ===========*)" << endl;
        if(!ramEP.next(kp, p, k3, x)) continue;
        subtract(k,kp,k1);
        set_v4(k2, ramEP.Pg);
        dbl_type Q2 = -mass2(k1);
        dbl_type s = sum_mass2(k1, k2);
        dbl_type t = subtract_mass2(k1, p);
        dbl_type u = subtract_mass2(k1, k3);

        println_4v_math("k", k);
        println_4v_math("P", P);
        cout << "(*--------*)" << endl;
        println_4v_math("k1", k1);
        println_4v_math("k2", k2);
        println_4v_math("kp", kp);
        cout << "(*--------*)" << endl;
        println_4v_math("p", p);
        println_4v_math("k3", k3);

        cout<<" Q2="<<Q2<<"; (* vs "<<ramEP.Q2<<"*)"<<endl;
        cout<<" s="<<s<<"; (* vs "<<sum_mass2(p,k3)<<"*)"<<endl;
        cout<<" t="<<t<<"; (* vs "<<subtract_mass2(k2, k3)<<"*)"<<endl;
        cout<<" u="<<u<<"; (* vs "<<subtract_mass2(k2, p)<<"*)"<<endl;

        assert(are_equal(Q2, ramEP.Q2));
        assert(are_equal(s, sum_mass2(p, k3)));
        assert(are_equal(t, subtract_mass2(k2, k3)));
        assert(are_equal(u, subtract_mass2(k2, p)));
        assert(are_equal(s + t + u, M2 - Q2));

        dbl_type alpha = 1. / 137, alphas = 0.3, PI = acos(-1);

        dbl_type matr2T = (-2097152 * alpha * pow(alphas, 2) * pow(M, 2) * pow(PI, 3) * pow(s, -2) * pow(Q2 + s, -2) * pow(u, -2)*
                ((s + u) * pow(Q2, 6) + (s + u) * pow(s, 4) * pow(u, 2) + u * pow(Q2, 4)*(-5 * s * u - 8 * pow(s, 2) + pow(u, 2)) + Q2 * u * pow(s, 3)*(4 * s * u + pow(s, 2) + 3 * pow(u, 2)) -
                pow(Q2, 5)*(5 * s * u + pow(s, 2) + 3 * pow(u, 2)) + s * pow(Q2, 3)*(-4 * u * pow(s, 2) + pow(s, 3) - s * pow(u, 2) + 3 * pow(u, 3)) +
                pow(Q2, 2) * pow(s, 2)*(u * pow(s, 2) + pow(s, 3) + 4 * s * pow(u, 2) + 4 * pow(u, 3))) * pow(s + u, -1) * pow(2 * Q2 + s + u, -2)) / 9.;
        cout << " matr2T=" << matr2T << endl;

        dbl_type matr2L = (-1048576 * alpha * Q2 * pow(alphas, 2) * pow(M, 2) * pow(PI, 3) * pow(s, -2) * pow(Q2 + s, -2) * pow(u, -2) * pow(s + u, -1)*
                (-2 * Q2 * u * (s + 2 * u) * pow(s, 3) - pow(s, 4) * pow(u, 2) - 2 * s * pow(Q2, 3)*(3 * s * u + 2 * pow(s, 2) + pow(u, 2)) - pow(Q2, 2) * pow(s, 2)*(8 * s * u + 2 * pow(s, 2) + 7 * pow(u, 2)) +
                2 * pow(Q2, 4) * pow(s + u, 2)) * pow(2 * Q2 + s + u, -2)) / 9.;
        cout << "matr2L=" << matr2L << endl;

        dbl_type matr2 = ramEP.wT * matr2T + ramEP.wL*matr2L;
        matr2=-matr2;
        cout << " matr2=" << matr2 << endl;
        assert(matr2>0);
    }
    return 0;
}
