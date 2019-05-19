#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_odeiv2.h>
#include "korb.h"

int main(void){
/*
This example program illustrates how you can incorporate the subroutines from 'korb.c' into your own programs

Appripriate commands for compiling:
--------------------------------

gcc korb.c -c -O2 -lm -lgsl -lgslcblas -lfftw3
gcc example.c -c -O2 -lm -lgsl -lgslcblas -lfftw3
gcc korb.o example.o -O2 -lm -lgsl -lgslcblas -lfftw3 -o example.out
rm -f korb.o
rm -f example.o
./example.out

--------------------------------

This should generate the executable 'example.out' that is capable of implementing any of the functions in 'korb.c'
*/
int eccentric = 1;
int inclined = 1;
double a = 0.5;
double e = 0.7;
double p = 7.0;
double x = 0.5;
double err = 1.0e-15;
korb_params orbpar;
// All orbital characteristics are calculated and stored inside "orbpar"
korb_getparams(eccentric,inclined,a,p,e,x,err,&orbpar);
// One useful application is to output the orbital characteristics
printf("E = %.15f\n", orbpar.E);
printf("Lz = %.15f\n", orbpar.Lz);
printf("Q = %.15f\n", orbpar.Q);
printf("Gamma = %.15f\n", orbpar.Ga);
printf("Omega_phi = %.15f\n", orbpar.wphi);
printf("Omega_th = %.15f\n", orbpar.wth);
printf("Omega_r = %.15f\n\n", orbpar.wr);
// Other useful applications include finding correlated quantities at certain times.
double lambda = 2.0;
// It is possible to treat Mino time as the independent variable
printf("t(lambda = %.15f) = %.15f\n", lambda, korb_tfromla(lambda,orbpar));
double chi_theta = M_PI/5.0;
// It is also possible to treat Mino time as the dependent variable
printf("lambda(chi_theta = %.15f) = %.15f\n\n", chi_theta, korb_lafromchi(chi_theta,orbpar));
// Because these quantities are constructed via Fourier series, it is inexpensive to repeat calculations at different times through the Fourier interpolant
lambda = 7.0;
chi_theta = 4*M_PI/5.0;
printf("t(lambda = %.15f) = %.15f\n", lambda, korb_tfromla(lambda,orbpar));
printf("lambda(chi_theta = %.15f) = %.15f\n\n", chi_theta, korb_lafromchi(chi_theta,orbpar));
// When you are finished extracting data associated with these orbital parameters, you should free any related memory to avoid memory leaks (especially if you are looping over a large number of orbits)
korb_freepar(orbpar);
// Once 'orbpar' is freed, you can reuse the variable for new cases
eccentric = 0;
inclined = 1;
a = 0.1;
e = 0.7;
p = 15.0;
x = -0.5;
err = 1.0e-15;
korb_getparams(eccentric,inclined,a,p,e,x,err,&orbpar);
printf("E = %.15f\n", orbpar.E);
printf("Lz = %.15f\n", orbpar.Lz);
printf("Q = %.15f\n", orbpar.Q);
printf("Gamma = %.15f\n", orbpar.Ga);
printf("Omega_phi = %.15f\n", orbpar.wphi);
printf("Omega_th = %.15f\n", orbpar.wth);
printf("Omega_r = %.15f\n", orbpar.wr);
korb_freepar(orbpar);
/* The original version of this example program has the following output:
-----------------------------------------------------------------

E = 0.966977955057975
Lz = 1.745837133883406
Q = 9.156020700594670
Gamma = 124.813457869459270
Omega_phi = 0.029819542436134
Omega_th = 0.027986782651302
Omega_r = 0.013422234455926

t(lambda = 2.000000000000000) = 301.838840775346057
lambda(chi_theta = 0.628318530717959) = 0.179906763530195

t(lambda = 7.000000000000000) = 920.257031562389557
lambda(chi_theta = 2.513274122871834) = 0.719457201206671


Because you made the 'eccentric' argument zero, the orbital eccentricity will be set to zero (regardless of your e input).

E = 0.969043201313842
Lz = -2.171102937060079
Q = 14.141521054487640
Gamma = 251.615327864292027
Omega_phi = -0.017197663429303
Omega_th = 0.017257492931553
Omega_r = 0.000000000000000

-----------------------------------------------------------------
*/
return 1;}


