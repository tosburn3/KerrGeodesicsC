// This software describes the generic bound motion of a test mass around a Kerr black hole. 
// Copyright (C) 2019  Thomas Osburn
//
// Please consider citing the following paper if you make use of this code:
// https://arxiv.org/abs/1905.13237
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

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

/*!
* Find z = cos^2(theta) from chi_theta (Eq. 2.26)
*/
double korb_zfromchi(double chi, korb_params orbpar){
double coschi = cos(chi);
return orbpar.zm*coschi*coschi;}

/*!
* Find theta from Mino time 
*/
double korb_thfromla(double la, korb_params orbpar){
if(orbpar.inclined==0)
  return M_PI/2;
double th = acos(sqrt(korb_zfromchi(korb_chifromla(la,orbpar),orbpar)));
if(orbpar.x<0)
  th = M_PI-th;
return th;}

/*!
* Find r from chi_r (Eq. 2.26)
*/
double korb_rfrompsi(double psi, korb_params orbpar){
return orbpar.p/(1+orbpar.e*cos(psi));}

/*!
* Find the derivative of Mino time with respect to chi_theta (Eq. 2.29)
*/
double korb_dladchi(double chi, korb_params orbpar){
return 1/sqrt(fabs(orbpar.B*(orbpar.zp-korb_zfromchi(chi,orbpar))));}

/*!
* Find the derivative of Mino time with respect to chi_r (Eq. 2.28)
*/
double korb_dladpsi(double psi, korb_params orbpar){
return (1-orbpar.e*orbpar.e)/sqrt(fabs((1-orbpar.E*orbpar.E)*(orbpar.p-orbpar.p3-orbpar.e*(orbpar.p+orbpar.p3*cos(psi)))*(orbpar.p-orbpar.p4+orbpar.e*(orbpar.p-orbpar.p4*cos(psi)))));}

/*!
* Find the derivative of chi_theta with respect to Mino time from chi_theta
*/
double korb_dchidla(double chi, korb_params orbpar){
return sqrt(fabs(orbpar.B*(orbpar.zp-korb_zfromchi(chi,orbpar))));}

/*!
* Find the derivative of chi_r with respect to Mino time from chi_r
*/
double korb_dpsidla(double psi, korb_params orbpar){
return sqrt(fabs((1-orbpar.E*orbpar.E)*(orbpar.p-orbpar.p3-orbpar.e*(orbpar.p+orbpar.p3*cos(psi)))*(orbpar.p-orbpar.p4+orbpar.e*(orbpar.p-orbpar.p4*cos(psi)))))/(1-orbpar.e*orbpar.e);}

/*! 
* Find Mino time from chi_r using a Fourier series (Eq. 3.5)
*/
double korb_lafrompsi(double psi, korb_params orbpar){
double la = orbpar.dladpsiamps[0]*psi;
int i;
for(i=1;i<orbpar.dladpsinum;i++)
  la += orbpar.dladpsiamps[i]/i*sin(i*psi);
return la;}

/*! 
* Find Mino time from chi_theta using a Fourier series
*/
double korb_lafromchi(double chi, korb_params orbpar){
double la = orbpar.dladchiamps[0]*chi;
int i;
for(i=1;i<orbpar.dladchinum;i++)
  la += orbpar.dladchiamps[i]/i*sin(i*chi);
return la;}

/*! 
* Find Delta from r
*/
double korb_D(double r, korb_params orbpar){
return r*r-2*r+orbpar.a*orbpar.a;}

/*! 
* Find T_r, which is the radial part of the derivative of t with respect to Mino time (Eq. 2.19)
*/
double korb_Tr(double psi, korb_params orbpar){
double r = korb_rfrompsi(psi,orbpar);
return orbpar.E*(r*r+orbpar.a*orbpar.a)*(r*r+orbpar.a*orbpar.a)/korb_D(r,orbpar)+orbpar.a*orbpar.Lz*(1-(r*r+orbpar.a*orbpar.a)/korb_D(r,orbpar));}

/*! 
* This function provides the integrand for calculating Mino time Fourier coefficients of T_r
*/
double complex korb_Trint(double psi, int n, korb_params orbpar){
double la = korb_lafrompsi(psi,orbpar);
return orbpar.Yr*korb_dladpsi(psi,orbpar)*korb_Tr(psi,orbpar)*(cos(n*orbpar.Yr*la)+I*sin(n*orbpar.Yr*la));}

/*! 
* Find T_theta, which is the polar part of the derivative of t with respect to Mino time (Eq. 2.20)
*/
double korb_Tth(double chi, korb_params orbpar){
return orbpar.E*orbpar.a*orbpar.a*(korb_zfromchi(chi,orbpar)-1);}

/*! 
* This function provides the integrand for calculating Mino time Fourier coefficients of T_theta
*/
double complex korb_Tthint(double chi, int n, korb_params orbpar){
double la = korb_lafromchi(chi,orbpar);
return orbpar.Yth*korb_dladchi(chi,orbpar)*korb_Tth(chi,orbpar)*(cos(n*orbpar.Yth*la)+I*sin(n*orbpar.Yth*la));}

/*! 
* Find Psi_r, which is the radial part of the derivative of phi with respect to Mino time (Eq. 2.17)
*/
double korb_Pr(double psi, korb_params orbpar){
double r = korb_rfrompsi(psi,orbpar);
return orbpar.a*orbpar.E*((r*r+orbpar.a*orbpar.a)/korb_D(r,orbpar)-1)-orbpar.a*orbpar.a*orbpar.Lz/korb_D(r,orbpar);}

/*! 
* This function provides the integrand for calculating Mino time Fourier coefficients of Psi_r
*/
double complex korb_Print(double psi, int n, korb_params orbpar){
double la = korb_lafrompsi(psi,orbpar);
return orbpar.Yr*korb_dladpsi(psi,orbpar)*korb_Pr(psi,orbpar)*(cos(n*orbpar.Yr*la)+I*sin(n*orbpar.Yr*la));}

/*! 
* Find Psi_theta, which is the polar part of the derivative of phi with respect to Mino time (Eq. 2.18)
*/
double korb_Pth(double chi, korb_params orbpar){
return orbpar.Lz/(1-korb_zfromchi(chi,orbpar));}

/*! 
* This function provides the integrand for calculating Mino time Fourier coefficients of Psi_theta
*/
double complex korb_Pthint(double chi, int n, korb_params orbpar){
double la = korb_lafromchi(chi,orbpar);
return orbpar.Yth*korb_dladchi(chi,orbpar)*korb_Pth(chi,orbpar)*(cos(n*orbpar.Yth*la)+I*sin(n*orbpar.Yth*la));}

/*! 
* This function provides the integrand for calculating Mino time Fourier coefficients of the derivative of chi_theta with respect to Mino time
*/
double complex korb_dchidlaint(double chi, int n, korb_params orbpar){
double la = korb_lafromchi(chi,orbpar);
return orbpar.Yth*(cos(n*orbpar.Yth*la)+I*sin(n*orbpar.Yth*la));}

/*! 
* This function provides the integrand for calculating Mino time Fourier coefficients of the derivative of chi_r with respect to Mino time
*/
double complex korb_dpsidlaint(double psi, int n, korb_params orbpar){
double la = korb_lafrompsi(psi,orbpar);
return orbpar.Yr*(cos(n*orbpar.Yr*la)+I*sin(n*orbpar.Yr*la));}

/*! 
* This function returns the integral of T_r with respect to Mino time given the Fourier coefficients of T_r (Eq. 3.19)
*/
double korb_dtrfromla(double la, korb_params orbpar){
double t = 0.0;
int j;
for(j=1;j<orbpar.Trnum;j++)
  t += orbpar.Tramps[j]/(j*orbpar.Yr)*sin(j*orbpar.Yr*la);
return t;}

/*! 
* This function returns the integral of T_theta with respect to Mino time given the Fourier coefficients of T_theta (Eq. 3.20)
*/
double korb_dtthfromla(double la, korb_params orbpar){
double t = 0.0;
int j;
for(j=1;j<orbpar.Tthnum;j++)
  t += orbpar.Tthamps[j]/(j*orbpar.Yth)*sin(j*orbpar.Yth*la);
return t;}

/*! 
* Find chi_r from Mino time
*/
double korb_psifromla(double la, korb_params orbpar){
if(orbpar.eccentric==0)
  return 0.0;  
double psi = orbpar.dpsidlaamps[0]*la;
int j;
for(j=1;j<orbpar.dpsidlanum;j++)
  psi += orbpar.dpsidlaamps[j]/(j*orbpar.Yr)*sin(j*orbpar.Yr*la);
return psi;}

/*! 
* Find chi_theta from Mino time
*/
double korb_chifromla(double la, korb_params orbpar){
if(orbpar.inclined==0)
  return 0.0;
double chi = orbpar.dchidlaamps[0]*la;
int j;
for(j=1;j<orbpar.dchidlanum;j++)
  chi += orbpar.dchidlaamps[j]/(j*orbpar.Yth)*sin(j*orbpar.Yth*la);
return chi;}

/*! 
* Find t from Mino time (Eq. 2.36)
*/
double korb_tfromla(double la, korb_params orbpar){
double t = orbpar.Ga*la;
if(orbpar.inclined!=0)
  t += korb_dtthfromla(la,orbpar);
if(orbpar.eccentric!=0)
  t += korb_dtrfromla(la,orbpar);
return t;}

/*! 
* This function returns the integral of Psi_r with respect to Mino time given the Fourier coefficients of Psi_r  (Eq. 3.21)
*/
double korb_dphirfromla(double la, korb_params orbpar){
double phi = 0.0;
int j;
for(j=1;j<orbpar.Prnum;j++)
  phi += orbpar.Pramps[j]/(j*orbpar.Yr)*sin(j*orbpar.Yr*la);
return phi;}

/*! 
* This function returns the integral of Psi_theta with respect to Mino time given the Fourier coefficients of Psi_theta  (Eq. 3.22)
*/
double korb_dphithfromla(double la, korb_params orbpar){
double phi = 0.0;
int j;
for(j=1;j<orbpar.Pthnum;j++)
  phi += orbpar.Pthamps[j]/(j*orbpar.Yth)*sin(j*orbpar.Yth*la);
return phi;}

/*! 
* Find phi from Mino time (Eq. 2.37)
*/
double korb_phifromla(double la, korb_params orbpar){
double phi = orbpar.Yphi*la;
if(orbpar.inclined!=0)
  phi += korb_dphithfromla(la,orbpar);
if(orbpar.eccentric!=0)
  phi += korb_dphirfromla(la,orbpar);
return phi;}

/*!
* Takes a function pointer and finds a single Fourier coefficient via DFT
*/
int korb_specint(double complex *amp, double large, int n, korb_params orbpar, double complex (*func)(double,int,korb_params)){
int N = 2;
double dx = 2*M_PI/N;
double complex sum = func(0,n,orbpar)+func(dx,n,orbpar);
*amp = sum/N;
double complex oldamp = 1.0e30;
double scale = fmax(large,cabs(*amp));
int i;
while(cabs(*amp-oldamp)/scale>orbpar.err||N<4*n)
{
oldamp = *amp;
for(i=0;i<N;i++)
  sum += func((i+.5)*dx,n,orbpar);
N *= 2;
dx = 2*M_PI/N;
*amp = sum/N;
scale = fmax(large,cabs(*amp));
}
return 1;}

/*!
* Takes a function pointer and finds a set of Fourier coefficients via DFT (including convergence assessment)
*/
int korb_dft(int *N, korb_params orbpar, double complex (*func)(double,int,korb_params), double **amps){
int num = 2;
double complex *camps;
camps = malloc(num * sizeof(double complex));
korb_specint(&(camps[0]),orbpar.err,0,orbpar,func);
double large = cabs(camps[0]);
korb_specint(&(camps[1]),large,1,orbpar,func);
camps[1] *= 2;
large = fmax(cabs(camps[0]),cabs(camps[1]));
double small = fmin(cabs(camps[0]),cabs(camps[1]));
while(small/large>orbpar.err||cabs(camps[num-1])/large>orbpar.err||cabs(camps[num-2])/large>orbpar.err||num<8)
{
num += 1;
camps = realloc(camps,num * sizeof(double complex));
korb_specint(&(camps[num-1]),large,num-1,orbpar,func);
camps[num-1] *= 2;
large = fmax(large,cabs(camps[num-1]));
small = fmin(small,cabs(camps[num-1]));
}
/* allow for non-even functions */ /*
*N = 2*(num-1);
*amps = malloc(*N * sizeof(double));
(*amps)[0] = creal(camps[0]);
(*amps)[*N/2] = creal(camps[num-1]);
int i;
for(i=1;i<*N/2;i++)
{
(*amps)[i] = creal(camps[i]);
(*amps)[*N-i] = cimag(camps[i]);
} */
/* even functions only */
*N = num;
*amps = malloc(*N * sizeof(double));
int i;
for(i=0;i<*N;i++)
  (*amps)[i] = creal(camps[i]);
/**/
free(camps);
return 1;}

/*!
* Takes a function pointer and samples it for DCT via FFTW
*/
int korb_dct(int num, double period, korb_params orbpar, double (*func)(double,korb_params), double amps[]){
fftw_plan p;
p = fftw_plan_r2r_1d(num,amps,amps,FFTW_REDFT00,FFTW_ESTIMATE);
double dx = period/(2*(num-1));
int n;
for(n=0;n<num;n++)
  amps[n] = func(n*dx,orbpar);
fftw_execute(p);
fftw_destroy_plan(p);
amps[0] *= 0.5;
amps[num-1] *= 0.5;
return 1;}

/*!
* Takes a function pointer and calculates DCT coefficients via FFTW (including convergence assessment)
*/
int korb_getamps(int *num, double period, korb_params orbpar, double (*func)(double,korb_params), double **amps){
*num = 4;
*amps = malloc(*num * sizeof(double));
int i;
for(i=0;i<*num;i++)
  (*amps)[i] = 1.0;
double scale = 1.0;
while(fabs((*amps)[*num-2])/scale>orbpar.err||fabs((*amps)[*num-1])/scale>orbpar.err)
  {
  *num = *num+2;
  *amps = realloc(*amps,*num * sizeof(double));
  korb_dct(*num,period,orbpar,func,*amps);
  scale = 0.0;
  for(i=0;i<*num;i++)
    scale = fmax(scale,fabs((*amps)[i]));
  }
double fact = 1.0/((*num)-1);
for(i=0;i<*num;i++)
  (*amps)[i] *= fact;
return 1;}

/*!
* Gives the derivative of the tortoise coordinate with respect to r (Eq. 2.50)
*/
double korb_drsdr(double rm, double a){
double rp = 1+sqrt(1-a*a);
return ((rm+rp)*(rm+rp)+a*a)/(rm*(rm+2*(rp-1)));}

/*!
* Gives the tortoise coordinate from r-r_+ (Eq. 2.49)
*/
double korb_rsfromrsubtrplus(double rsubtrplus, double a){
double rp = 1+sqrt(1-a*a);
return rsubtrplus+rp+(rp*log(.5*rsubtrplus)+(rp-2)*log(.5*rsubtrplus+rp-1))/(rp-1);}

/*!
* Residual between guessing the tortoise coordinate value and the actual value (for root finding)
*/
double korb_rmf(double rm, void *params){
double a = ((double *)params)[0];
double rs = ((double *)params)[1];
return korb_rsfromrsubtrplus(rm,a)-rs;}

/*!
* Gives the derivative of the tortoise coordinate with respect to r (for root finding)
*/
double korb_rmdf(double rm, void *params){
double a = ((double *)params)[0];
return  korb_drsdr(rm,a);}

/*!
* Gives both the residual between guessing the tortoise coordinate value and the actual value and the derivative of the tortoise coordinate with respect to r (for root finding)
*/
void korb_rmfdf(double rm, void *params, double *y, double *dy){
double a = ((double *)params)[0];
double rs = ((double *)params)[1];
*y = korb_rsfromrsubtrplus(rm,a)-rs;
*dy = korb_drsdr(rm,a);}

/*!
* Uses various initial guess strategies to efficiently and accurately find r-r_+ from the tortoise coordinate via root finding
*/
double korb_rsubtrplusfromrs(double rs, double a){
double ars[2] = {a,rs};
double rp = 1+sqrt(1-a*a);
double guess;
if(rs<-2)
  guess = 2*pow(1-a*a,1/rp-.5)*exp((a*a-rp+(rp-1)*rs)/rp);
else if(rs<1000)
  guess = 2*gsl_sf_lambert_W0(exp(.5*(rs-rp)));
else
  guess = rs;
int status;
int iter = 0;
int max_iter = 1000;
const gsl_root_fdfsolver_type *T;
gsl_root_fdfsolver *s;
double x0, x = guess;
gsl_function_fdf FDF;
FDF.f = &korb_rmf;
FDF.df = &korb_rmdf;
FDF.fdf = &korb_rmfdf;
FDF.params = ars;
T = gsl_root_fdfsolver_steffenson;
s = gsl_root_fdfsolver_alloc(T);
gsl_root_fdfsolver_set(s,&FDF,x);
do
  {
  iter++;
  status = gsl_root_fdfsolver_iterate(s);
  x0 = x;
  x = gsl_root_fdfsolver_root(s);
  status = gsl_root_test_delta(x,x0,fmax(1.0,rs)*3.0e-16,3.0e-16);
  }
while(status==GSL_CONTINUE&&iter<max_iter);
gsl_root_fdfsolver_free(s);
return x;}

/*!
* This function gathers all necessary numerical techniques to determine generic orbital trajectories of a test mass around a Kerr black hole. The general strategy involves spectral calculation of Fourier coefficients describing the rate-of-change of position, then the position itself is determined by integrating the Fourier series term-by-term. All necessary information for reconstructing the position at a certain time is encoded in the "korb_params" data structure. It adaptively handles convergence testing and memory allocation (but a separate function is needed to free the allocated memory).
*/ 
int korb_getparams(int eccentric, int inclined, double a, double p, double e0, double x0, double err, korb_params *orbpar){
double e, x;
// Checks are performed to ensure consistency between the boolean flags indicating the presence of eccentricity and/or inclination with the floating point values of eccentricity and/or inclination. Other checks are also made to exclude scenarios that are unphysical or beyond the capabilities of this code.
if(eccentric==0)
{
if(fabs(e0)>1.0e-14)
  {
  printf("\nBecause you made the 'eccentric' argument zero, the orbital eccentricity will be set to zero (regardless of your e input).\n\n");
  }
orbpar->eccentric = 0;
e = 5.0e-16;
}
else
{
if(fabs(e0)<1.0e-14)
  {
  printf("\nBecause you input such a small (near zero) orbital eccentricity value, the 'eccentric' parameter will be set to zero.\n\n");
  orbpar->eccentric = 0;
  e = 5.0e-16;
  }
else if(fabs(e0)>0.9+5.0e-16)
  {
  printf("\nWarning, eccentricities larger than 0.9 not supported, setting e = 0.9\n\n");
  orbpar->eccentric = 1;
  e = 0.9;
  }
else if(e0<0)
  {
  printf("\nWarning, negative eccentricity, the eccentricity will be made positive.\n\n");
  orbpar->eccentric = 1;
  e = fabs(e0);
  }
else
  {
  orbpar->eccentric = 1;
  e = e0;
  }
}
if(inclined==0)
{
if(fabs(1-fabs(x0))>1.0e-14)
  {
  printf("\nBecause you made the 'inclined' argument zero, this orbit is assumed to be equatorial (regardless of your x input). Prograde or retrograde implied from the sign of your x argument.\n\n");
  }
if(x0>=0)
  x = cos(5.0e-16);
else
  x = cos(M_PI-5.0e-16);
orbpar->inclined = 0;
}
else
{
if(fabs(1-fabs(x0))<=1.0e-14)
  {
  printf("\nBecause you input such a small (near zero) inclination angle value (|x| is near 1), the 'inclined' parameter will be set to zero. Prograde or retrograde implied from the sign of your x argument.\n\n");
  if(x0>0)
    x = cos(1.0e-15);
  else
    x = cos(M_PI-5.0e-16);
  orbpar->inclined = 0;
  }
else if(x0>1+1.0e-14)
  {
  printf("\nWarning, x > 1 is not physical, resetting x = 1 and setting the 'inclined' parameter to zero.\n\n");
  x = cos(1.0e-15);
  orbpar->inclined = 0;
  }
else if(x0<-1-1.0e-14)
  {
  printf("\nWarning, x < -1 is not physical, resetting x = -1 and setting the 'inclined' parameter to zero.\n\n");
  x = cos(M_PI-5.0e-16);
  orbpar->inclined = 0;
  }
else if(fabs(x0)<0.15)
  {
  printf("\nWarning, near polar orbits are not supported. Setting |x| = 0.15, prograde or retrograde implied from the sign of your x argument.\n\n");
  if(x0>0)
    x = 0.15;
  else
    x = -0.15;
  orbpar->inclined = 1;
  }
else
  {
  orbpar->inclined = 1;
  x = x0;
  }
}
if(err<1.0e-15)
  orbpar->err = 1.0e-15;
else
  orbpar->err = err;
orbpar->a = a;
orbpar->e = e;
orbpar->p = p;
orbpar->x = x;
orbpar->r1 = p/(1-e);
orbpar->r2 = p/(1+e);
double zm = 1-x*x;
double r1 = orbpar->r1;
double r2 = orbpar->r2;
// The mathematical relationships governing the constants of motion are hard-coded
orbpar->E = sqrt(fabs((2*a*a*(1-e*e)*(1-e*e)*x*x*(p*p-a*a*(1-e*e)*(1-x*x))*(-(p*p*p*(-4+3*p+e*e*(4+p)))+a*a*(-1+e*e)*(-1+e*e)*p*p*(-2+x*x)+a*a*a*a*(-1+e*e)*(-1+e*e)*(-1+e*e)*(-1+x*x))-2*a*(-1+e*e)*(-1+e*e)*x*(p*p+a*a*(1-e*e)*(-1+x*x))*sqrt(fabs(p*(a*a*a*a*(-1+e*e)*(-1+e*e)+(-4*e*e+(-2+p)*(-2+p))*p*p+2*a*a*p*(-2+p+e*e*(2+p)))*(p*p+a*a*(1-e*e)*(-1+x*x))*(p*p*p*p-2*a*a*(1+e*e)*p*p*(-1+x*x)+a*a*a*a*(-1+e*e)*(-1+e*e)*(-1+x*x)*(-1+x*x))))+((-4*e*e+(-2+p)*(-2+p))*p*p*p+a*a*a*a*(-1+e*e)*(-1+e*e)*(-1+x*x)*(-p+(-1+e*e+p)*x*x)+a*a*p*p*(2*(-2+p+e*e*(2+p))-(-3+e*e*e*e+2*p+2*e*e*(1+p))*x*x))*(p*p*p*p*(-3-e*e+p)+a*a*a*a*(-1+e*e)*(-1+e*e)*(-1+x*x)*(-1+e*e-p+(-1+e*e+p)*x*x)-2*a*a*p*p*(1+e*e*e*e+p*(-1+x*x)+e*e*(-2+p*(-1+x*x)))))/(4*a*a*(1-e*e)*(1-e*e)*x*x*(p*p-a*a*(1-e*e)*(1-x*x))*((-3-e*e)*p*p*p*p+a*a*(-1+e*e)*(-1+e*e)*p*p*(-2+x*x)+a*a*a*a*(-1+e*e)*(-1+e*e)*(-1+e*e)*(-1+x*x))+pow((3+e*e-p)*p*p*p*p-a*a*a*a*(-1+e*e)*(-1+e*e)*(-1+x*x)*(-1+e*e-p+(-1+e*e+p)*x*x)+2*a*a*p*p*(1+e*e*e*e+p*(-1+x*x)+e*e*(-2+p*(-1+x*x))),2))));
orbpar->Lz = (2*a*(1+e)*(1+e)*p*x*x*orbpar->E-x*sqrt(fabs((a*a*(1+e)*(1+e)+p*(-2-2*e+p))*(-p*p+a*a*(1+e)*(1+e)*(-1+x*x))*(a*a*(1+e)*(1+e)*(-1+x*x)*(-1+orbpar->E*orbpar->E)+p*(-2-2*e+p-p*orbpar->E*orbpar->E)))))/((1+e)*((2+2*e-p)*p+a*a*(1+e)*(1+e)*(-1+x*x)));
orbpar->Q = (p*((1+e)*(1+e)*orbpar->Lz*orbpar->Lz*(2+2*e-p)-4*a*(1+e)*(1+e)*(1+e)*orbpar->Lz*orbpar->E+a*a*(1+e)*(1+e)*(-p+(2+2*e+p)*orbpar->E*orbpar->E)+p*p*(2+2*e+p*(-1+orbpar->E*orbpar->E))))/(a*a*(1+e)*(1+e)*(1+e)*(1+e)+(1+e)*(1+e)*p*(-2-2*e+p));
orbpar->B = a*a*(1-orbpar->E*orbpar->E);
orbpar->zm = 1-x*x;
orbpar->thmin = acos(sqrt(fabs(orbpar->zm)));
orbpar->zp = (orbpar->B+orbpar->Lz*orbpar->Lz+orbpar->Q+sqrt(fabs(-4*orbpar->B*orbpar->Q+pow(orbpar->B+orbpar->Lz*orbpar->Lz+orbpar->Q,2))))/(2*orbpar->B);
orbpar->r3 = (2*a*a-orbpar->B*(r1+r2)+sqrt(fabs(r1*r2*(-4*a*a*a*a*orbpar->B*orbpar->Q+r1*r2*pow(-2*a*a+orbpar->B*(r1+r2),2))))/(r1*r2))/(2*orbpar->B);
orbpar->r4 = -(-2*a*a*r1*r2+orbpar->B*r1*r1*r2+orbpar->B*r1*r2*r2+sqrt(fabs(r1*r2*(-4*a*a*a*a*orbpar->B*orbpar->Q+r1*r2*pow(-2*a*a+orbpar->B*(r1+r2),2)))))/(2*orbpar->B*r1*r2);
orbpar->p3 = orbpar->r3*(1-e);
orbpar->p4 = orbpar->r4*(1+e);
orbpar->Yphi = 0.0;
orbpar->Ga = 0.0;
// When appropriate, Fourier amplitudes governing orbital evolution are calculated adaptively
if(orbpar->inclined!=0)
  {
  korb_getamps(&(orbpar->dladchinum),2*M_PI,*orbpar,korb_dladchi,&(orbpar->dladchiamps));
  orbpar->Vth = 2*M_PI*orbpar->dladchiamps[0];
  orbpar->Yth = 2*M_PI/orbpar->Vth;
  korb_dft(&(orbpar->Tthnum),*orbpar,korb_Tthint,&(orbpar->Tthamps));
  korb_dft(&(orbpar->Pthnum),*orbpar,korb_Pthint,&(orbpar->Pthamps));
  orbpar->Yphi += orbpar->Pthamps[0];
  orbpar->Ga += orbpar->Tthamps[0];
  korb_dft(&(orbpar->dchidlanum),*orbpar,korb_dchidlaint,&(orbpar->dchidlaamps));
  }
if(orbpar->eccentric!=0)
  {
  korb_getamps(&(orbpar->dladpsinum),2*M_PI,*orbpar,korb_dladpsi,&(orbpar->dladpsiamps));
  orbpar->Vr = 2*M_PI*orbpar->dladpsiamps[0];
  orbpar->Yr = 2*M_PI/orbpar->Vr;
  korb_dft(&(orbpar->Trnum),*orbpar,korb_Trint,&(orbpar->Tramps));
  korb_dft(&(orbpar->Prnum),*orbpar,korb_Print,&(orbpar->Pramps));
  orbpar->Yphi += orbpar->Pramps[0];
  orbpar->Ga += orbpar->Tramps[0];
  korb_dft(&(orbpar->dpsidlanum),*orbpar,korb_dpsidlaint,&(orbpar->dpsidlaamps));
  }
if(orbpar->inclined==0)
  {
  orbpar->Yphi += korb_Pth(0.0,*orbpar);
  orbpar->Ga += korb_Tth(0.0,*orbpar);
  orbpar->wth = 0;
  }
if(orbpar->eccentric==0)
  {
  orbpar->Yphi += korb_Pr(0.0,*orbpar);
  orbpar->Ga += korb_Tr(0.0,*orbpar);
  orbpar->wr = 0;
  }
orbpar->wphi = orbpar->Yphi/orbpar->Ga;
if(orbpar->inclined!=0)
  orbpar->wth = orbpar->Yth/orbpar->Ga;
if(orbpar->eccentric!=0)
  orbpar->wr = orbpar->Yr/orbpar->Ga;
return 1;}

/*!
* This function frees the memory allocated by the function "korb_getparams"
*/ 
int korb_freepar(korb_params orbpar){
if(orbpar.eccentric==1)
{
free(orbpar.dladpsiamps);
free(orbpar.dpsidlaamps);
free(orbpar.Tramps);
free(orbpar.Pramps);
}
if(orbpar.inclined==1)
{
free(orbpar.dladchiamps);
free(orbpar.dchidlaamps);
free(orbpar.Tthamps);
free(orbpar.Pthamps);
}
return 1;}

/*!
* The main function takes arguments: eccentric (boolean), inclined (boolean), a (double), p (double), e (double), x (double), lambdaSteps (integer), lambdaMax (double), and from those arguments calculates and prints the constants of motion, orbital frequencies, and position at equally spaced Mino times. Must compile with flag '-DLINK_DEFAULT' when not linking through another program.
*/ 
#ifdef LINK_DEFAULT
int main(int argc, char **argv){
int eccentric;
int inclined;
double a;
double e;
double p;
double x;
int num;
double lambdamax;
// While the code is capable of sacrificing accuracy in favor of a bit more efficiency, in this main function the relative error is fixed at a minimal value (near the limitations of double precision)
double err = 1.0e-15;
// Parse inputs and check for conflicts
if(argc!=9)
  {
  printf("Wrong number of inputs.\n\n");
  printf("The 8 necessary inputs are:\n");
  printf(" 1. 'eccentric' (boolean)\n");
  printf(" 2. 'inclined' (boolean)\n");
  printf(" 3. 'a' (decimal)\n");
  printf(" 4. 'p' (decimal)\n");
  printf(" 5. 'e' (decimal)\n");
  printf(" 6. 'x' (decimal)\n");
  printf(" 7. 'lambdaNum' (integer)\n");
  printf(" 8. 'lambdaMax' (decimal)\n\n");
  return 0;
  }
eccentric = strtol(argv[1],NULL,0);
inclined = strtol(argv[2],NULL,0);
a = strtod(argv[3],NULL);
p = strtod(argv[4],NULL);
e = strtod(argv[5],NULL);
x = strtod(argv[6],NULL);
num = strtol(argv[7],NULL,0);
if(num<1)
{
num = 1;
printf("lambdaSteps must be at least 1, setting lambdaSteps = 1\n");
}
lambdamax = strtod(argv[8],NULL);
if(lambdamax<0)
{
lambdamax = 10.0;
printf("lambdaMax must be positive, setting lambdaMax = 10.0\n");
}
korb_params orbpar;
// All orbital characteristics are calculated and stored inside "orbpar"
korb_getparams(eccentric,inclined,a,p,e,x,err,&orbpar);
// Relevant orbital characteristics within "orbpar" are printed
printf("E = %.15f\n", orbpar.E);
printf("Lz = %.15f\n", orbpar.Lz);
printf("Q = %.15f\n", orbpar.Q);
printf("Gamma = %.15f\n", orbpar.Ga);
printf("Omega_phi = %.15f\n", orbpar.wphi);
printf("Omega_th = %.15f\n", orbpar.wth);
printf("Omega_r = %.15f\n", orbpar.wr);
// Equipped with the information inside "orbpar", positions are reconstructed and printed at various Mino time values based on the user inputs
int j;
double delta_lambda = lambdamax/num;
for(j=1;j<=num;j++)
{
printf("\n\n   lambda = %.15f\n", j*delta_lambda);
printf("t = %.15f\n", korb_tfromla(j*delta_lambda,orbpar));
printf("r = %.15f\n", korb_rfrompsi(korb_psifromla(j*delta_lambda,orbpar),orbpar)); 
printf("theta = %.15f\n", korb_thfromla(j*delta_lambda,orbpar));
printf("phi = %.15f\n", korb_phifromla(j*delta_lambda,orbpar));
printf("chi_r = %.15f\n", korb_psifromla(j*delta_lambda,orbpar));
printf("chi_theta = %.15f\n", korb_chifromla(j*delta_lambda,orbpar));
}
// Finally, the memory allocated for "orbpar" is freed
korb_freepar(orbpar);
return 1;}
#endif


