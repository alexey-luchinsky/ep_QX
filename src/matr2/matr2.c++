#include "kinematics/RamboEP.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>

extern dbl_type mc, Mcc, ecm, PI, alpha, alphas, x;
extern dbl_type Opsi, O3S11, O3S18, O3P08, O1S08;
extern dbl_type P[4], k[4], kp[4], k1[4], k2[4], k3[4], pPsi[4];
extern bool gauge;
extern bool check;
cmplx epsG1[4], epsG2[4], epsG3[4], epsPsi[4];
extern dbl_type Q2;

void write_histogram_to_file(TH1F &histogram, string file_name) {
    const char *__file_name__ = file_name.c_str();
    remove(__file_name__);
    ofstream file;

    file.open(__file_name__);
    for (int i = 1; i <= histogram.GetNbinsX(); i++)
        file << setiosflags(ios::scientific) << histogram.GetBinCenter(i) <<
        " " << setiosflags(ios::scientific) << histogram.GetBinContent(i) / histogram.GetBinWidth(i) <<
        " " << setiosflags(ios::scientific) << histogram.GetBinError(i) / histogram.GetBinWidth(i) << endl;

    file.close();
}

void print_debug() {
    cout<<" ecm="<<ecm<<";"<<endl;
    println_4v_math("P",P);
    println_4v_math("k",k);
    println_4v_math("kp",kp);
    println_4v_math("k3",k3);
    println_4v_math("pPsi",pPsi);
    cout<<" (kp+pPsi+k3)^2="<<sum_mass2(kp, pPsi, k3)<<" vs "<<x*ecm*ecm<<endl;
    assert(are_equal(sum_mass2(kp, pPsi, k3), x * ecm * ecm));
};

bool kinematics(RamboEP *ramEP) {
    bool ok=true;
    if (!ramEP->next(kp, pPsi, k3, x)) return false;
    set_v4(k, ramEP->kIn); // k
    ok = ok && is_zero(mass2(k));
//    assert(is_zero(mass2(k)));
    set_v4(kp, ramEP->kOut); // kp
    ok = ok && is_zero(mass2(kp));
//    assert(is_zero(mass2(kp)));
    set_v4(k2, ramEP->Pg); // k2
    ok = ok && are_equal(sum_mass2(pPsi, k3), ramEP->W2);
//    assert(are_equal(sum_mass2(pPsi, k3), ramEP->W2));
    subtract(k, kp, k1); // k1=k-kp
    ok = ok && are_equal(mass2(k1), -ramEP->Q2);
    ok = ok && are_equal(mass2(pPsi), Mcc * Mcc);
    ok = ok && is_zero(mass2(k3));
    ok = ok && are_equal(sum_mass2(k, k2), x * ecm * ecm);
    ok = ok && are_equal(sum_mass2(kp, pPsi, k3), x * ecm * ecm);
//    assert(are_equal(mass2(k1), -ramEP->Q2));
//    assert(are_equal(mass2(pPsi), Mcc * Mcc));
//    assert(is_zero(mass2(k3)));
//    assert(are_equal(sum_mass2(k, k2), x * ecm * ecm));
//    assert(are_equal(sum_mass2(kp, pPsi, k3), x * ecm * ecm));
    if(!ok && check) {
        print_debug();
        return false;
    };
    return true;
};

cmplx amp_3S1_cs() {
    cmplx sunDelta=COMPLEX_ZERO;
    sunDelta += sp(epsG1,epsPsi)*((sp(epsG2,k3)) * ((sp(epsG3,pPsi)) * (sp(k2,pPsi)) + -1 * ((sp(epsG3,k2)) * (sp(k2,k3) + sp(k2,pPsi) + -1 * (sp(k3,pPsi))))) + (sp(epsG2,epsG3)) * ((sp(k2,k3) + sp(k2,pPsi)) * (sp(k2,k3) + -1 * (sp(k3,pPsi)))) + (sp(epsG2,pPsi)) * (-1 * ((sp(epsG3,pPsi)) * (sp(k2,k3))) + (sp(epsG3,k2)) * (sp(k3,pPsi))));
    sunDelta += sp(epsG1,epsG3)*(-1 * ((sp(k2,k3)) * (-1 * ((sp(epsG2,pPsi)) * (sp(epsPsi,k3))) + (sp(epsG2,epsPsi)) * (sp(k2,k3) + sp(k2,pPsi)))) + (sp(epsG2,k3)) * (-1 * ((sp(epsPsi,k3)) * (sp(k2,pPsi))) + (sp(epsPsi,k2)) * (sp(k2,k3) + sp(k2,pPsi))));
    sunDelta += sp(epsG1,epsG2)*(-1 * ((sp(epsG3,pPsi)) * ((sp(epsPsi,k2)) * (sp(k2,k3)))) + (sp(epsG3,epsPsi)) * ((sp(k2,k3)) * (-1 * (sp(k2,k3)) + sp(k3,pPsi))) + (sp(epsG3,k2)) * ((sp(epsPsi,k3)) * (sp(k2,k3) + -1 * (sp(k3,pPsi))) + (sp(epsPsi,k2)) * (sp(k3,pPsi))));
    sunDelta += sp(epsG1,k3)*(-1 * ((sp(epsG2,pPsi)) * ((sp(epsG3,epsPsi)) * (sp(k2,k3)))) + -1 * ((sp(epsG2,epsG3)) * ((sp(epsPsi,k2)) * (sp(k2,k3)))) + (sp(epsG2,k3)) * ((sp(epsG3,epsPsi)) * (sp(k2,pPsi))) + -1 * ((sp(epsG2,epsG3)) * ((sp(epsPsi,k2)) * (sp(k2,pPsi)))) + (sp(epsG2,epsPsi)) * ((sp(epsG3,k2)) * (sp(k2,k3) + sp(k2,pPsi))));
    sunDelta += sp(epsG1,k2)*((sp(epsG2,k3)) * ((sp(epsG3,epsPsi)) * (sp(k2,k3) + -1 * (sp(k3,pPsi)))) + (sp(epsG2,epsG3)) * ((sp(epsPsi,k3)) * (-1 * (sp(k2,k3)) + sp(k3,pPsi))) + (sp(epsG2,epsPsi)) * ((sp(epsG3,pPsi)) * (sp(k2,k3)) + -1 * ((sp(epsG3,k2)) * (sp(k3,pPsi)))));
    sunDelta *= alphas * (alpha * (-128/9. * ((pow(Q2,-1)) * ((pow(PI,2)) * ((pow(sp(k2,pPsi),-1)) * ((pow(2 * (sp(k2,k3)) + sp(k2,pPsi) + -1 * (sp(k3,pPsi)),-1)) * ((pow(sp(k3,pPsi),-1)) * (sqrt(mc * O3S11)))))))));
    sunDelta *= sqrt(8.); // deta_ab^2 
    return sunDelta;
}

cmplx amp_3S1_co() {
	 cmplx sund=COMPLEX_ZERO;
	 sund += sp(epsG1,epsPsi)*((sp(epsG2,k3)) * ((sp(epsG3,pPsi)) * (sp(k2,pPsi)) + -1 * ((sp(epsG3,k2)) * (sp(k2,k3) + sp(k2,pPsi) + -1 * (sp(k3,pPsi))))) + (sp(epsG2,epsG3)) * ((sp(k2,k3) + sp(k2,pPsi)) * (sp(k2,k3) + -1 * (sp(k3,pPsi)))) + (sp(epsG2,pPsi)) * (-1 * ((sp(epsG3,pPsi)) * (sp(k2,k3))) + (sp(epsG3,k2)) * (sp(k3,pPsi))));
	 sund += sp(epsG1,epsG3)*(-1 * ((sp(k2,k3)) * (-1 * ((sp(epsG2,pPsi)) * (sp(epsPsi,k3))) + (sp(epsG2,epsPsi)) * (sp(k2,k3) + sp(k2,pPsi)))) + (sp(epsG2,k3)) * (-1 * ((sp(epsPsi,k3)) * (sp(k2,pPsi))) + (sp(epsPsi,k2)) * (sp(k2,k3) + sp(k2,pPsi))));
	 sund += sp(epsG1,epsG2)*(-1 * ((sp(epsG3,pPsi)) * ((sp(epsPsi,k2)) * (sp(k2,k3)))) + (sp(epsG3,epsPsi)) * ((sp(k2,k3)) * (-1 * (sp(k2,k3)) + sp(k3,pPsi))) + (sp(epsG3,k2)) * ((sp(epsPsi,k3)) * (sp(k2,k3) + -1 * (sp(k3,pPsi))) + (sp(epsPsi,k2)) * (sp(k3,pPsi))));
	 sund += sp(epsG1,k3)*(-1 * ((sp(epsG2,pPsi)) * ((sp(epsG3,epsPsi)) * (sp(k2,k3)))) + -1 * ((sp(epsG2,epsG3)) * ((sp(epsPsi,k2)) * (sp(k2,k3)))) + (sp(epsG2,k3)) * ((sp(epsG3,epsPsi)) * (sp(k2,pPsi))) + -1 * ((sp(epsG2,epsG3)) * ((sp(epsPsi,k2)) * (sp(k2,pPsi)))) + (sp(epsG2,epsPsi)) * ((sp(epsG3,k2)) * (sp(k2,k3) + sp(k2,pPsi))));
	 sund += sp(epsG1,k2)*((sp(epsG2,k3)) * ((sp(epsG3,epsPsi)) * (sp(k2,k3) + -1 * (sp(k3,pPsi)))) + (sp(epsG2,epsG3)) * ((sp(epsPsi,k3)) * (-1 * (sp(k2,k3)) + sp(k3,pPsi))) + (sp(epsG2,epsPsi)) * ((sp(epsG3,pPsi)) * (sp(k2,k3)) + -1 * ((sp(epsG3,k2)) * (sp(k3,pPsi)))));
	 sund *= alphas * (alpha * (-32/3. * ((pow(Q2,-1)) * ((pow(PI,2)) * ((pow(sp(k2,pPsi),-1)) * ((pow(2 * (sp(k2,k3)) + sp(k2,pPsi) + -1 * (sp(k3,pPsi)),-1)) * ((pow(sp(k3,pPsi),-1)) * ((sqrt(2/3.)) * (sqrt(mc * O3S18))))))))));
	 sund *= sqrt(40/3.); // d_abc^2 
         return sund;
}

cmplx amp_1S0_co() {
    cmplx sunf=COMPLEX_ZERO;
    sunf += lcv(epsG1,epsG3,k2,k3)*((sp(epsG2,pPsi)) * ((sp(k2,k3)) * ((-1) * (sp(k2,k3)) + sp(k3,pPsi))));
    sunf += lcv(epsG1,epsG3,k3,pPsi)*((-1) * (((sp(epsG2,pPsi)) * (sp(k2,k3)) + (-1) * ((sp(epsG2,k3)) * (sp(k2,pPsi)))) * (sp(k2,k3) + (-1) * (sp(k3,pPsi)))));
    sunf += lcv(epsG1,epsG2,k2,pPsi)*((-1) * ((sp(k2,k3) + sp(k2,pPsi)) * ((sp(epsG3,pPsi)) * (sp(k2,k3)) + (-1) * ((sp(epsG3,k2)) * (sp(k3,pPsi))))));
    sunf += lcv(epsG1,epsG2,epsG3,pPsi)*(pow(sp(k2,k3),3) + (sp(k2,k3)) * ((sp(k2,pPsi)) * (sp(k3,pPsi))));
    sunf += lcv(epsG1,epsG2,k2,k3)*((sp(epsG3,pPsi)) * ((sp(k2,pPsi)) * ((-1) * (sp(k2,k3)) + sp(k3,pPsi))));
    sunf += lcv(epsG1,epsG2,epsG3,k3)*((sp(k2,k3)) * ((sp(k2,pPsi)) * ((-1) * (sp(k2,k3)) + sp(k3,pPsi))));
    sunf += lcv(epsG1,epsG2,epsG3,k2)*((-1) * ((sp(k2,k3)) * ((sp(k2,k3) + sp(k2,pPsi)) * (sp(k3,pPsi)))));
    sunf += sp(epsG1,epsG3)*((lcv(epsG2,k2,k3,pPsi)) * (pow(sp(k2,k3),2) + (sp(k2,pPsi)) * (sp(k3,pPsi))));
    sunf += lcv(epsG1,epsG3,k2,pPsi)*((sp(epsG2,k3)) * (pow(sp(k2,k3),2) + (sp(k2,pPsi)) * (sp(k3,pPsi))));
    sunf *= (-32/3.) * ((alpha) * ((alphas) * ((pow(mc,-1/2.)) * ((pow(Q2,-1)) * ((pow(PI,2)) * ((pow(sp(k2,k3),-1)) * ((pow(sp(k2,pPsi),-1)) * ((pow((2) * (sp(k2,k3)) + sp(k2,pPsi) + (-1) * (sp(k3,pPsi)),-1)) * ((pow(sp(k3,pPsi),-1)) * (sqrt(O1S08)))))))))));
    sunf *= sqrt(24); // f_abc^2
    return sunf;
}

cmplx amp_3P0_co() {
    cmplx sunf=COMPLEX_ZERO;
    sunf += sp(epsG1,epsG2)*(-1 * ((sp(k2,pPsi)) * (((sp(epsG3,pPsi)) * (sp(k2,k3)) + -1 * ((sp(epsG3,k2)) * (sp(k3,pPsi)))) * (8 * ((pow(mc,2)) * (pow(sp(k2,k3),3))) + (sp(k2,k3)) * ((3 * (pow(sp(k2,pPsi),2)) + pow(sp(k3,pPsi),2) + 4 * ((sp(k2,pPsi)) * (3 * (pow(mc,2)) + -1 * (sp(k3,pPsi))))) * (sp(k3,pPsi))) + (sp(k2,pPsi)) * ((sp(k2,pPsi) + -1 * (sp(k3,pPsi))) * ((4 * (pow(mc,2)) + sp(k2,pPsi) + -1 * (sp(k3,pPsi))) * (sp(k3,pPsi)))) + (pow(sp(k2,k3),2)) * (-1 * (pow(sp(k3,pPsi),2)) + (sp(k2,pPsi)) * (4 * (pow(mc,2)) + 3 * (sp(k3,pPsi))))))));
    sunf += sp(epsG1,epsG3)*(-1 * (((sp(epsG2,pPsi)) * (sp(k2,k3)) + -1 * ((sp(epsG2,k3)) * (sp(k2,pPsi)))) * ((sp(k3,pPsi)) * (8 * ((pow(mc,2)) * (pow(sp(k2,k3),3))) + (sp(k2,pPsi)) * ((sp(k2,pPsi) + -1 * (sp(k3,pPsi))) * ((4 * (pow(mc,2)) + sp(k2,pPsi) + -1 * (sp(k3,pPsi))) * (sp(k3,pPsi)))) + -1 * ((pow(sp(k2,k3),2)) * (pow(sp(k2,pPsi),2) + 4 * ((pow(mc,2)) * (sp(k3,pPsi))) + -3 * ((sp(k2,pPsi)) * (sp(k3,pPsi))))) + -1 * ((sp(k2,k3)) * ((sp(k2,pPsi)) * (pow(sp(k2,pPsi),2) + -4 * ((sp(k2,pPsi)) * (sp(k3,pPsi))) + 3 * ((sp(k3,pPsi)) * (-4 * (pow(mc,2)) + sp(k3,pPsi))))))))));
    sunf += sp(epsG1,pPsi)*((sp(k2,k3)) * ((sp(k2,pPsi)) * ((sp(k3,pPsi)) * ((sp(k2,pPsi) + sp(k3,pPsi)) * ((sp(epsG2,k3)) * (-1 * ((sp(epsG3,pPsi)) * (sp(k2,pPsi))) + (sp(epsG3,k2)) * (sp(k2,k3) + sp(k2,pPsi) + -1 * (sp(k3,pPsi)))) + -1 * ((sp(epsG2,epsG3)) * ((sp(k2,k3) + sp(k2,pPsi)) * (sp(k2,k3) + -1 * (sp(k3,pPsi))))) + (sp(epsG2,pPsi)) * ((sp(epsG3,pPsi)) * (sp(k2,k3)) + -1 * ((sp(epsG3,k2)) * (sp(k3,pPsi)))))))));
    sunf += sp(epsG1,k3)*((sp(epsG2,pPsi)) * ((sp(epsG3,pPsi)) * ((sp(k2,k3)) * ((sp(k2,pPsi)) * ((pow(sp(k3,pPsi),2)) * (-4 * (pow(mc,2)) + -1 * (sp(k2,pPsi)) + sp(k3,pPsi)) + (pow(sp(k2,k3),2)) * (-8 * (pow(mc,2)) + 2 * (sp(k3,pPsi))) + -2 * ((sp(k2,k3)) * ((sp(k2,pPsi)) * (2 * (pow(mc,2)) + -1 * (sp(k3,pPsi))) + (sp(k3,pPsi)) * (-4 * (pow(mc,2)) + sp(k3,pPsi))))))) + (sp(epsG3,k2)) * ((sp(k3,pPsi)) * (8 * ((pow(mc,2)) * (pow(sp(k2,k3),3))) + (pow(sp(k2,pPsi),2)) * ((4 * (pow(mc,2)) + sp(k2,pPsi) + -1 * (sp(k3,pPsi))) * (sp(k3,pPsi))) + -1 * ((sp(k2,k3)) * ((sp(k2,pPsi)) * (pow(sp(k2,pPsi),2) + (sp(k3,pPsi)) * (-4 * (pow(mc,2)) + sp(k3,pPsi)) + -2 * ((sp(k2,pPsi)) * (2 * (pow(mc,2)) + sp(k3,pPsi)))))) + (pow(sp(k2,k3),2)) * (-1 * (pow(sp(k2,pPsi),2)) + -4 * ((pow(mc,2)) * (sp(k3,pPsi))) + (sp(k2,pPsi)) * (8 * (pow(mc,2)) + sp(k3,pPsi)))))) + (sp(k2,pPsi)) * (-1 * ((sp(epsG2,epsG3)) * ((sp(k2,k3) + sp(k2,pPsi)) * ((sp(k3,pPsi)) * (8 * ((pow(mc,2)) * (pow(sp(k2,k3),2))) + (sp(k2,pPsi)) * ((4 * (pow(mc,2)) + sp(k2,pPsi) + -1 * (sp(k3,pPsi))) * (sp(k3,pPsi))) + (sp(k2,k3)) * ((sp(k2,pPsi)) * (4 * (pow(mc,2)) + -1 * (sp(k2,pPsi)) + sp(k3,pPsi))))))) + (sp(epsG2,k3)) * (4 * ((pow(mc,2)) * ((sp(epsG3,k2)) * ((sp(k2,k3)) * ((sp(k3,pPsi)) * (sp(k2,pPsi) + sp(k3,pPsi)))))) + (sp(epsG3,pPsi)) * ((sp(k2,pPsi)) * ((pow(sp(k2,k3),2)) * (8 * (pow(mc,2)) + -2 * (sp(k3,pPsi))) + (pow(sp(k3,pPsi),2)) * (4 * (pow(mc,2)) + sp(k2,pPsi) + -1 * (sp(k3,pPsi))) + 2 * ((sp(k2,k3)) * ((sp(k2,pPsi)) * (2 * (pow(mc,2)) + -1 * (sp(k3,pPsi))) + (sp(k3,pPsi)) * (-4 * (pow(mc,2)) + sp(k3,pPsi)))))))));
    sunf += sp(epsG1,k2)*(-1 * ((sp(k3,pPsi)) * ((sp(epsG2,epsG3)) * ((sp(k2,pPsi)) * ((sp(k2,k3) + -1 * (sp(k3,pPsi))) * (8 * ((pow(mc,2)) * (pow(sp(k2,k3),2))) + (sp(k2,pPsi)) * ((4 * (pow(mc,2)) + sp(k2,pPsi) + -1 * (sp(k3,pPsi))) * (sp(k3,pPsi))) + -1 * ((sp(k2,k3)) * ((sp(k3,pPsi)) * (4 * (pow(mc,2)) + -1 * (sp(k2,pPsi)) + sp(k3,pPsi))))))) + -1 * ((sp(epsG2,pPsi)) * (((sp(epsG3,pPsi)) * (sp(k2,k3)) + -1 * ((sp(epsG3,k2)) * (sp(k3,pPsi)))) * (2 * ((pow(sp(k2,k3),2)) * (4 * (pow(mc,2)) + sp(k2,pPsi))) + (pow(sp(k2,pPsi),2)) * (4 * (pow(mc,2)) + sp(k2,pPsi) + -1 * (sp(k3,pPsi))) + 2 * ((sp(k2,k3)) * (pow(sp(k2,pPsi),2) + (sp(k2,pPsi)) * (4 * (pow(mc,2)) + -1 * (sp(k3,pPsi))) + -2 * ((pow(mc,2)) * (sp(k3,pPsi)))))))))) + (sp(epsG2,k3)) * ((sp(k2,pPsi)) * (-4 * ((pow(mc,2)) * ((sp(epsG3,k2)) * ((sp(k2,k3)) * ((sp(k3,pPsi)) * (sp(k2,pPsi) + sp(k3,pPsi)))))) + (sp(epsG3,pPsi)) * (8 * ((pow(mc,2)) * (pow(sp(k2,k3),3))) + (pow(sp(k3,pPsi),2)) * ((sp(k2,pPsi)) * (-4 * (pow(mc,2)) + -1 * (sp(k2,pPsi)) + sp(k3,pPsi))) + (sp(k2,k3)) * ((sp(k3,pPsi)) * (pow(sp(k2,pPsi),2) + (sp(k2,pPsi)) * (4 * (pow(mc,2)) + -2 * (sp(k3,pPsi))) + (sp(k3,pPsi)) * (4 * (pow(mc,2)) + sp(k3,pPsi)))) + (pow(sp(k2,k3),2)) * ((sp(k2,pPsi)) * (4 * (pow(mc,2)) + sp(k3,pPsi)) + -1 * ((sp(k3,pPsi)) * (8 * (pow(mc,2)) + sp(k3,pPsi))))))));
    sunf *= alphas * (alpha * (-32/3. * ((pow(3,-1/2.)) * ((pow(mc,-3/2.)) * ((pow(Q2,-1)) * ((pow(PI,2)) * ((pow(sp(k2,k3),-1)) * ((pow(sp(k2,pPsi),-2)) * ((pow(2 * (sp(k2,k3)) + sp(k2,pPsi) + -1 * (sp(k3,pPsi)),-2)) * ((pow(sp(k3,pPsi),-2)) * (sqrt(O3P08))))))))))));
    sunf *= sqrt(24); // f_abc^2
    return sunf;
}

dbl_type getMatr2(cmplx (*ampFunc)(), int nPolPsi) {
    dbl_type matr2 = 0;
    dbl_type Q2 = -mass2(k1);
    for (int iG1 = 0; iG1 < 2; ++iG1) {
        lepton_current(iG1, k, kp, epsG1);
        for (int iG2 = 0; iG2 < 2; ++iG2) {
            set_gluon_polarization(iG2, k2, epsG2);
            if (gauge) for (int i = 0; i < 4; ++i) epsG2[i] = k2[i];
            for (int iG3 = 0; iG3 < 2; ++iG3) {
                set_gluon_polarization(iG3, k3, epsG3);
                for (int iPsi = 0; iPsi < nPolPsi; ++iPsi) {
                    set_psi_polarization(iPsi, pPsi, epsPsi);
                    cmplx amp = (*ampFunc)();
                    matr2 += pow(abs(amp), 2);
                };
            };
        };
    };
    return matr2;
}

dbl_type getMatr2_3S1_cs() {
    return getMatr2(&amp_3S1_cs, 3);
}

dbl_type getMatr2_3S1_co() {
    return getMatr2(&amp_3S1_co, 3);
}

dbl_type getMatr2_1S0_co() {
    return getMatr2(&amp_1S0_co, 1);
}

dbl_type getMatr2_3P0_co() {
    return getMatr2(&amp_3P0_co, 1);
}



