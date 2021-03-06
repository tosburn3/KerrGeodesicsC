# Kerr Geodesics in C

This software describes the generic bound motion of a test mass around a Kerr black hole. There are also functions for calculating r from the tortoise coordinate via root finding.

The original implementation of this code was in support of the following work:
https://arxiv.org/abs/1905.13237

Please consider citing the above paper if you make use of this code. All descriptions involving numbered equations in this code's documentation are in reference to the above paper.

### Dependencies

KerrGeodesicsC depends on:

- The GNU Scientific Library (version 2, https://www.gnu.org/software/gsl/)
- FFTW (version 3, http://www.fftw.org/)

### Compiling the code

The code can typically be compiled in gcc with the following command:

gcc korb.c -O2 -lm -lgsl -lgslcblas -lfftw3 -o korb.out -DLINK_DEFAULT

The flag '-DLINK_DEFAULT' should be omitted when linking through another program with its own 'main' function. Ideally this will generate the executable 'korb.out' for calculating Kerr orbital trajectories. However, you may have to add additional flags in the event that any dependent libraries or headers are located in non-standard directories.

It is also possible to use this code as a library of functions that can be compiled and linked with other programs, see 'example.c' for some simple demsonstrations.

### Usage

Geodesics can be computed with the command:

./korb.out eccentric inclined a p e x lambdaNum lambdaMax

Consider the following example:

./korb.out 1 0 0.99 7.0 0.9 0.5 2 1.0

These parameters describe an eccentric (eccentric = 1) equatorial (inclined = 0) orbit with a/M = 0.99, p = 7.0, e = 0.9, that is prograde (because x > 0, the specific x value is ignored when inclined = 0). Position values will be printed at 2 Mino times (lambdaNum = 2) up to a maximum Mino time of 1.0 (lambdaMax = 1.0). The output generated from these parameters is:

______________________________________________
Because you made the 'inclined' argument zero, this orbit is assumed to be equatorial (regardless of your x input). Prograde or retrograde implied from the sign of your x argument.

E = 0.986565156231870
Lz = 3.052859990816571
Q = 0.000000000000004
Gamma = 568.585782351074272
Omega_phi = 0.005924590732463
Omega_th = 0.000000000000000
Omega_r = 0.004150454187797


   lambda = 0.500000000000000
t = 14.391447729962465
r = 4.940154465458621
theta = 1.570796326794897
phi = 1.792448368140982
chi_r = 1.089093839765892
chi_theta = 0.000000000000000


   lambda = 1.000000000000000
t = 67.984126677545021
r = 17.159459542406832
theta = 1.570796326794897
phi = 3.453882116286580
chi_r = 2.288752031370152
chi_theta = 0.000000000000000
____________________________________________

### Licence

The code is licensed under the GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)




