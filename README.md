# ep_QX

In this package I try to calculate $J/\psi$ meson photoproduction 
$$e(k) p(P) \to e(k') J/\psi(p_\psi) +X$$
at different colliders (D0, LHep, etc).

The package consists of the following parts

## kinematics/

Useful subroutines for kinematics and phase space generation.

__Random.h, Random.cpp__

Random generator

__algebra.h, algebra.cpp__

Some functions working with 4-momenta, etc

__Rambo2.h, Rambo2.cpp__

Class for 2D lorentz-invariant phase space (LIPS) generation. I use [Okun] agreement for $n$-particle LIPS:
$$
d\Phi_n[P;p_1,\dots p_n] = (2\pi)^4\delta^4\left(P-\sum_i p_i\right) \prod_i \frac{d^3p_i}{(2\pi)^3 2E_i}.
$$
In the case of two massless particles it correspond to integrated LIPS
$$
\Phi_2 = \frac{1}{8\pi}
$$

The main methods of the class are:

```Rambo2(dbl_type m1, dbl_type m2, Random *random)```: constructor

``` dbl_type next(dbl_type ecm, dbl_type (&k1)[4], dbl_type (&k2)[4])```: event generation. Returns weight of the event and fills the momenta of the outgoing particles. These momenta satisfy the relation
$$ k_1 + k_2 = \{\text{ecm};0,0,0\}$$

Thus, the code
```c++
#include "kinematics/Rambo2.h"
int main(void) {
	Random *random=new Random();
	Rambo2 rambo(0,0,random);
	dbl_type ecm=10, k1[4], k2[4];
	dbl_type sum=0;
	int nEv=1e4;
	for(int iEv=0; iEv<nEv; ++iEv)
		sum += rambo(ecm,k1,k2)/nEv;
	cout<<sum<<endl;
}
```
should return something like
$$ \frac{1}{8\pi} \approx 0.0398$$