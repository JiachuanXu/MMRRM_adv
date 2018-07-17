/*
	FLAGS, Constants and PATH
*/
#ifndef _MMRRM_ADV_
#define _MMRRM_ADV_

//includes standard h
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <complex.h>
#include <unistd.h>
#include <fftw3.h>
#ifdef _OPENMP
	#include <omp.h>
#endif

#ifdef _OPENMP
	#define nthreads (int)4// Number of threads
	#define mem 8		   // Total memory allocated
#else
	#define nthreads (int)1
	#define mem 4
#endif

/*	Define constants	*/
#define PI (double)(3.14159265358979323846264338327)
#define E  (double)(2.71828182846)
#define YEAR_START 1900
#define MON_START 1
/*	Define auxiliary macros	*/
#define box_vol(dim) (unsigned long long)((dim)*(dim)*(dim))	//The volume of the box, except for velocity boxes
#define kbox_vol(dim) (unsigned long long)((dim)*(dim)*((dim)/2+1llu)) 	//The volume of the FFT transformed boxes
#define in_tr(i,j,k,dim) (unsigned long long)(((dim)*(i)+(j))*(dim)+(k))	//Index transformation for boxes. e.g. in_tr(1,1,1)=dim^2+dim+1
#define in_tr_c(i,j,k,dim) (unsigned long long)(k+(dim/2+1llu)*(j+dim*i)) //Index transformation for complex box
#define bin_index(kint,imu) (unsigned long long)((kint + ORDER/2 + 1)*(kint - ORDER/2)/2 + imu) // Index transformation for mu decomposition

/****  Begin user changable definitions ****/
#define M_FLOAT  35	//The max lenth for a float number in character format
#define M_PATH  300	//The max lenth for boxes_PATH
#define M_BOXNAME  300	//The max lenth for box name
#define M_CMD 600	//The max lenth for command line
#define M_TEMP 300
#define ARG_C 11 //Number of terms is argument
#define k_factor 1.4	//k_factor for PS
#define kmincut (int)6 // kmincut should >= 2
#define ORDER (int)4// The highest order for mu decomposition fitting
#define std_arg "../prmt/args_mmrrm_adv.dat"	//The standard dir of arguments

/*	Define cosmological constant and parameters	*/
#define a(z) (double)(1.0/(1.0+(z)))
#define h (double)(0.678)
#define H (double)(100 * (h))	//Hubble constant H0, in km.s^-1.Mpc^-1
#define Om  (double)(0.308)	//Omega matter
#define Ol  (double)(1.0-Om)	//Omega lambda
#define Ob  (double)(0.048)	//Omega baryons
#define Ez(z) (double)(sqrt(Om * pow((double)(1.0+z),3.0) + Ol)) // Matter dominant
#define Hz(z) (double)(H * Ez(z))// Matter dominant
#define cell_size(box_size, dim) (double)((box_size)/(dim))	//In units of Mpc
#define T21 (double)(0.068) //temperature corresponding to the 21cm photon, K
#define T_cmb0  (double)(2.728) // K
#define T_cmb(z)  T_cmb0*(1.0+(double)z) // K
#define A10  (double)(2.85e-15) // spontaneous emission coefficient in s^-1
#define C  (double)(29979245800.0) // speed of light cm/s
#define NU0  (double)(1420406000.0) //Frequency of 21cm emission, Hz, s^-1
#define Mpc_over_km  (double)(3.086e19) //1 Mpc = NU0 km
#define CRE_TAU  (double)(0.1)  // Cretical optical depth
#define G_const (double)(6.67408e-11)	//in standard unit, m^3*kg^-1*s^-2
#define rho_cre (double)(3*H*H/(8*PI*G_const*1000*Mpc_over_km*Mpc_over_km))	//Cretical density today, in g*cm^-3
#define Y_He (double)(0.245)	//Helium fraction
#define MU_H (double)(1.67353276e-24)	//Mass of one H atom, g
#define Ts_lim (double)(9.999e10) //Set the Ts in high temperature limit

/****	FLAGS	****/
/* Flag for observer frame of output data cube  */
#define MESH2MESH 0 // 1 to calculate observables in redshift space, 0 in real space
/* Flag for approximations in RSD correction */
/* If you do not have Ts information, set both OPTHIN & HIGHTS = 1*/
#define OPTHIN 1 // 1 for artificially applying optical thin approximation, 0 for standard routine
#define HIGHTS 1 // 1 for artificially applying Ts>>T_cmb approximation, 0 for standard routine
/* Flags for analysis subroutines */
// You could just turn off all of them if you only need the output data cube
//#define QLIN 0 // 1 to turn on quasi-linear \mu_k decomposition, 0 to turn off. (see Y. Mao er al. 2012)
#define POWINK 1 // 1 to turn on calculating power spectrum in k interval, 0 for not
#define POWINLOGK 1 // 1 to turn on calculating power spectrum in log(k) interval, 0 for not
#define MUDECINK 1 // 1 to turn on mu decomposition in k interval, 0 to turn off
#define MUDECINLOGK 1 // 1 to turn of mu decomposition in log(k) interval, 0 to turn off
#define TBPDF 1 // 1 to turn on the pdf calculation of Tb, 0 to turn off
#define STATIS 1	// 1 to turn on the calculation of statistics of sknewness, kurtoisis, mean, dev, rms
#define GLOBAL_EVOL 1 // 1 to record the global evolution
#define DEBUG 0 // 1 to turn on debugging 

/*	PATH 	*/
#define MMRRM_BOX_OP "../Output_boxes"// Output dir of MMRRM boxes
#define MMRRM_PS_OP "../Output_PS" // Output dir of MMRRM power spectrum and other statistics


#endif
