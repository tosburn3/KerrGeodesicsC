/*! \mainpage  korb.h
 *  \brief     This software describes the generic bound motion of a test mass around a Kerr black hole. There are also functions for calculating r from the tortoise coordinate via root finding. The main function takes inputs: eccentric (boolean), inclined (boolean), a (double), p (double), e (double), x (double), lambdaSteps (integer), lambdaMax (double), and from those arguments calculates and prints the constants of motion, orbital frequencies, and position at equally spaced Mino times.
 *  \author    Thomas Osburn
 *  \version   0.1
 *  \date      2019
 *  \copyright 2019 Thomas Osburn
 */

/*
Copyright 2019 Thomas Osburn
*/

/*!
* This data structure encodes all information needed to reconstruct the position of a test mass orbiting a Kerr black hole at any time via Fourier series. The pointers within will become lists of Fourier coefficients once they are allocated by the function "korb_getparams"
*/ 
typedef struct korb_params korb_params;
struct korb_params{
int eccentric; /*!< Boolean to indicate whether the orbit is eccentric */
int inclined; /*!< Boolean to indicate whether the orbit is inclined */
double a; /*!< specific spin angular momentum divided by black hole mass */
double p; /*!< Orbital semilatus rectum */
double e; /*!< Orbital eccentricity */
double x; /*!< x = sin(theta_min) */
double thmin; /*!< Minimum polar angle */
double E; /*!< Specific energy */
double Lz; /*!< Specific angular momentum */
double Q; /*!< Carter constant */
double B; /*!< Beta = (1-E^2)*a^2 */
double zm;  /*!< z_- = cos^2(theta_min) */
double zp; /*!< z_+ = other root of polar equation of motion (not directly related to theta_max) */
double r1; /*!< r_1 = largest root of radial equation of motion = r_max */
double r2; /*!< r_2 = 2nd largest root of radial equation of motion = r_min */
double r3; /*!< r_3 = 3rd largest root of radial equation of motion (not directly related to a turning point) */
double r4; /*!< r_4 = smallest root of radial equation of motion (not directly related to a turning point) */
double p3; /*!< p_3 = r_3*(1-e) */
double p4; /*!< p_4 = r_4*(1+e) */
double Vr; /*!< Lambda_r = radial period in Mino time */
double Vth; /*!< Lambda_theta = polar period in Mino time */
double Yr; /*!< Upsilon_r = fundamental radial angular frequency in Mino time */
double Yth; /*!< Upsilon_theta = fundamental polar angular frequency in Mino time */
double Yphi; /*!< Upsilon_phi = fundamental azimuthal angular frequency in Mino time */
double Ga; /*!< Gamma = average rate of t advance in Mino time */
double wr; /*!< Omega_r = fundamental radial angular frequency in t */
double wth; /*!< Omega_theta = fundamental polar angular frequency in t */
double wphi; /*!< Omega_phi = fundamental azimuthal angular frequency in t */
double *dladpsiamps;/*!< Fourier amplitudes of the derivative of Mino time with respect to chi_r */
int dladpsinum; /*!< number of Fourier amplitudes of the derivative of Mino time with respect to chi_r */
double *dladchiamps; /*!< Fourier amplitudes of the derivative of Mino time with respect to chi_theta */
int dladchinum; /*!< number of Fourier amplitudes of the derivative of Mino time with respect to chi_theta */
double *dpsidlaamps; /*!< Fourier amplitudes of the derivative of chi_r with respect to Mino time */
int dpsidlanum; /*!< number of Fourier amplitudes of the derivative of chi_r with respect to Mino time */
double *dchidlaamps; /*!< Fourier amplitudes of the derivative of chi_theta with respect to Mino time */
int dchidlanum; /*!< number of Fourier amplitudes of the derivative of chi_theta with respect to Mino time */
double *Tramps; /*!< Fourier amplitudes of T_r */
int Trnum; /*!< number of Fourier amplitudes of T_r */
double *Tthamps; /*!< Fourier amplitudes of T_theta */
int Tthnum; /*!< number of Fourier amplitudes of T_theta */
double *Pramps; /*!< Fourier amplitudes of Psi_r */
int Prnum; /*!< number of Fourier amplitudes of Psi_r */
double *Pthamps; /*!< Fourier amplitudes of Psi_theta */
int Pthnum; /*!< number of Fourier amplitudes of Psi_theta */
double err; /*!< relative error tolerance used to determine convergence of each Fourier series */
};

/*!
* This function gathers all necessary numerical techniques to determine generic orbital trajectories of a test mass around a Kerr black hole. The general strategy involves spectral calculation of Fourier coefficients describing the rate-of-change of position, then the position itself is determined by integrating the Fourier series term-by-term. All necessary information for reconstructing the position at a certain time is encoded in the "korb_params" data structure. It adaptively handles convergence testing and memory allocation (but a separate function is needed to free the allocated memory).
*/
extern int korb_getparams(int eccentric, int inclined, double a, double p, double e0, double x0, double err, korb_params *orbpar);

/*!
* This function frees the memory allocated by the function "korb_getparams"
*/
extern int korb_freepar(korb_params orbpar);

/*!
* Find z = cos^2(theta) from chi_theta
*/
extern double korb_zfromchi(double chi, korb_params orbpar);

/*!
* Find theta from Mino time
*/
extern double korb_thfromla(double la, korb_params orbpar);

/*!
* Find r from chi_r
*/
extern double korb_rfrompsi(double psi, korb_params orbpar);

/*!
* Find the derivative of Mino time with respect to chi_theta
*/
extern double korb_dladchi(double chi, korb_params orbpar);

/*!
* Find the derivative of Mino time with respect to chi_r
*/
extern double korb_dladpsi(double psi, korb_params orbpar);

/*!
* Find the derivative of chi_theta with respect to Mino time from chi_theta
*/
extern double korb_dchidla(double chi, korb_params orbpar);

/*!
* Find the derivative of chi_r with respect to Mino time from chi_r
*/
extern double korb_dpsidla(double psi, korb_params orbpar);

/*!
* Find the derivative of chi_r with respect to Mino time from chi_r
*/
extern double korb_lafrompsi(double psi, korb_params orbpar);

/*! 
* Find Mino time from chi_theta using a Fourier series
*/
extern double korb_lafromchi(double chi, korb_params orbpar);

/*! 
* Find Delta from r
*/
extern double korb_D(double r, korb_params orbpar);

/*! 
* Find T_r, which is the radial part of the derivative of t with respect to Mino time
*/
extern double korb_Tr(double psi, korb_params orbpar);

/*! 
* This function provides the integrand for calculating Mino time Fourier coefficients of T_r
*/
extern double complex korb_Trint(double psi, int n, korb_params orbpar);

/*! 
* Find T_th, which is the polar part of the derivative of t with respect to Mino time
*/
extern double korb_Tth(double chi, korb_params orbpar);

/*! 
* This function provides the integrand for calculating Mino time Fourier coefficients of T_th
*/
extern double complex korb_Tthint(double chi, int n, korb_params orbpar);

/*! 
* Find Psi_r, which is the radial part of the derivative of phi with respect to Mino time
*/
extern double korb_Pr(double psi, korb_params orbpar);

/*! 
* This function provides the integrand for calculating Mino time Fourier coefficients of Psi_r
*/
extern double complex korb_Print(double psi, int n, korb_params orbpar);

/*! 
* Find Psi_theta, which is the polar part of the derivative of phi with respect to Mino time
*/
extern double korb_Pth(double chi, korb_params orbpar);

/*! 
* This function provides the integrand for calculating Mino time Fourier coefficients of Psi_theta
*/
extern double complex korb_Pthint(double chi, int n, korb_params orbpar);

/*! 
* This function provides the integrand for calculating Mino time Fourier coefficients of the derivative of chi_theta with respect to Mino time
*/
extern double complex korb_dchidlaint(double chi, int n, korb_params orbpar);

/*! 
* This function provides the integrand for calculating Mino time Fourier coefficients of the derivative of chi_r with respect to Mino time
*/
extern double complex korb_dpsidlaint(double psi, int n, korb_params orbpar);

/*! 
* This function returns the integral of T_r with respect to Mino time given the Fourier coefficients of T_r ( called orbpar.Tramps[] )
*/
extern double korb_dtrfromla(double la, korb_params orbpar);

/*! 
* This function returns the integral of T_theta with respect to Mino time given the Fourier coefficients of T_theta ( called orbpar.Tthamps[] )
*/
extern double korb_dtthfromla(double la, korb_params orbpar);

/*! 
* Find chi_r from Mino time
*/
extern double korb_psifromla(double la, korb_params orbpar);

/*! 
* Find chi_theta from Mino time
*/
extern double korb_chifromla(double la, korb_params orbpar);

/*! 
* Find t from Mino time
*/
extern double korb_tfromla(double la, korb_params orbpar);

/*! 
* This function returns the integral of Psi_r with respect to Mino time given the Fourier coefficients of Psi_r ( called orbpar.Pramps[] )
*/
extern double korb_dphirfromla(double la, korb_params orbpar);

/*! 
* This function returns the integral of Psi_theta with respect to Mino time given the Fourier coefficients of Psi_theta ( called orbpar.Pthamps[] )
*/
extern double korb_dphithfromla(double la, korb_params orbpar);

/*! 
* Find phi from Mino time
*/
extern double korb_phifromla(double la, korb_params orbpar);

/*!
* Takes a function pointer and finds a single Fourier coefficient via DFT
*/
extern int korb_specint(double complex *amp, double large, int n, korb_params orbpar, double complex (*func)(double,int,korb_params));

/*!
* Takes a function pointer and finds a set of Fourier coefficients via DFT (including convergence assessment)
*/
extern int korb_dft(int *N, korb_params orbpar, double complex (*func)(double,int,korb_params), double **amps);

/*!
* Takes a function pointer and samples it for DCT via FFTW
*/
extern int korb_dct(int num, double period, korb_params orbpar, double (*func)(double,korb_params), double amps[]);

/*!
* Takes a function pointer and calculates DCT coefficients via FFTW (including convergence assessment)
*/
extern int korb_getamps(int *num, double period, korb_params orbpar, double (*func)(double,korb_params), double **amps);

/*!
* Gives the derivative of the tortoise coordinate with respect to r
*/
extern double korb_drsdr(double rm, double a);

/*!
* Gives the tortoise coordinate from r-r_+
*/
extern double korb_rsfromrsubtrplus(double rsubtrplus, double a);

/*!
* Residual between guessing the tortoise coordinate value and the actual value (for root finding)
*/
extern double korb_rmf(double rm, void *params);

/*!
* Gives the derivative of the tortoise coordinate with respect to r (for root finding)
*/
extern double korb_rmdf(double rm, void *params);

/*!
* Gives both the residual between guessing the tortoise coordinate value and the actual value and the derivative of the tortoise coordinate with respect to r (for root finding)
*/
extern void korb_rmfdf(double rm, void *params, double *y, double *dy);

/*!
* Uses various initial guess strategies to efficiently and accurately find r-r_+ from the tortoise coordinate via root finding
*/
extern double korb_rsubtrplusfromrs(double rs, double a);

