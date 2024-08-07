/* voigt.c: computes exact Voigt profiles in terms of the complex
   error function (requires libcerf) */

#include "mex.h"
#include <cerf.h>
#include <math.h>
#include <stdio.h>

#define LAMBDAS_ARG          prhs[0]
#define Z_ARG                prhs[1]
#define N_ARG                prhs[2]
#define NUM_LINES_ARG        prhs[3]
#define sigma_ARG            prhs[4]
#define pixel_sigma_ARG      prhs[5]
#define PROFILE_ARG          plhs[0]

/* number of lines in CIV doublet */
#define NUM_LINES 2

/* note: all units are CGS */

/* physical constants */

static const double c   = 2.99792458e+10;              /* speed of light          cm s⁻¹        */
// static const double k   = 1.38064852e-16;         /* Boltzmann constant      erg K⁻¹       */
// static const double m_p = 1.672621898e-24;        /* proton mass             g             */
// static const double m_e = 9.10938356e-28;         /* electron mass           g             */
// double e = 1.6021766208e-19 * c / 10; */
// static const double e   = 4.803204672997660e-10;  /* elementary charge       statC         */

/* CIV doublet */

static const double transition_wavelengths[] =         /* transition wavelengths  cm            */
  {
    1.5481949462890625e-5,
    1.55077001953125e-5   
      };

static const double oscillator_strengths[] =           /* oscillator strengths    dimensionless */
  {
    0.189900,
    0.094750
  };

static const double Gammas[] =                         /* transition rates        s^-1          */
  {
    2.643e+08,
    2.628e+08
  };

/* assumed constant */
/* static const double T = 1e+04; */                   /* gas temperature         K             */

/* derived constants */

/* b = sqrt(2 * k * T / m_p); */
/* b = sqrt(2 * k * T / (6*m_p+6*m_n)); */
/* static const double b =
/*     1.28486551932562422e+06; */                       /* Doppler parameter       cm s⁻¹        */

/* sigma = b / M_SQRT2; */
// static const double sigma = 262181.7413311349;         /* T=1e4 */
// static const double sigma= 370780.9743970856;            /*  T=2e4 */
//static const double sigma= 829091.4635154926;            /* T=1e5*/
// static const double sigma= 2621817.413311349;            /* T=1e6*/
// leading_constants[i] =
      // M_PI * e * e * oscillator_strengths[i] * transition_wavelengths[i] / (m_e * c) ;
// */
static const double leading_constants[] =              /* leading constants  cm² */
  {     7.802895118381213e-08,
     3.899701297867750e-08
 };

/* gammas[i] = Gammas[i] * transition_wavelengths[i] / (4 * M_PI); */
static const double gammas[] =                         /* Lorentzian widths       cm s⁻¹        */
  {
     3.255002952981575e+02,
     3.243136695286643e+02
  };

static const int width = 3;                      /* width of convolution     dimensionless */

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]) {

  double *lambdas;
  double *profile;
  double *multipliers;
  double *raw_profile;
  double *instrument_profile;
  double *pixel_sigma;
  double z;
  double N;
  double velocity;
  double total;
  double sigma;
  int num_lines, i, j, k, ii, jj;
  mwSize num_points;

  /* get input */
  lambdas     = mxGetPr(LAMBDAS_ARG);                /* wavelengths             Å             */
  z           = mxGetScalar(Z_ARG);                  /* redshift                dimensionless */
  N           = mxGetScalar(N_ARG);                  /* column density          cm⁻²          */
  sigma       = mxGetScalar(sigma_ARG);              /* Doppler broadening  paramberter   cm/s*/
  pixel_sigma = mxGetPr(pixel_sigma_ARG);             /* pixel sigma from spSpec                */

  num_lines = (nrhs > 3) ? (int)(mxGetScalar(NUM_LINES_ARG)) : NUM_LINES;
  num_points = mxGetNumberOfElements(LAMBDAS_ARG);
  /* initialize output */
  PROFILE_ARG = mxCreateDoubleMatrix(num_points - 2 * width, 1, mxREAL);
  profile = mxGetPr(PROFILE_ARG);                /* absorption profile      dimensionless */

  /* to hold the profile before instrumental broadening */
  raw_profile = mxMalloc(num_points * sizeof(double));
  multipliers = mxMalloc(num_lines * sizeof(double));
  instrument_profile = mxMalloc((2*width+1)* sizeof(double));

  for (i = 0; i < num_lines; i++)
    multipliers[i] = c / (transition_wavelengths[i] * (1 + z)) / 1e8;

  /* compute raw Voigt profile */
  for (i = 0; i < num_points; i++) {
    /* apply each absorption line */
    total = 0;
    for (j = 0; j < num_lines; j++) {
      /* velocity relative to transition wavelength */
      velocity = lambdas[i] * multipliers[j] - c;
      total += -leading_constants[j] * voigt(velocity, sigma, gammas[j]);
    }

    raw_profile[i] = exp(N * total);
  }

  num_points = mxGetNumberOfElements(PROFILE_ARG);

  /* instrumental broadening */
  for (i = 0; i < num_points; i++){
    if ((i>=width) && (i<=num_points-width)){
      total = 0;
      for (ii = -width, jj = 0; ii <= width; ii++, jj++) {
        
        instrument_profile[jj] = exp(-0.5 * ii * ii / (pixel_sigma[i] * pixel_sigma[i]))/pixel_sigma[i];
        total += instrument_profile[jj];
        }
      
      for (ii = 0; ii < 2 * width + 1; ii++)
        instrument_profile[ii] /= total;

      for (j = i, k = 0; j <= i + 2 * width; j++, k++){
        profile[i] += raw_profile[j] * instrument_profile[k];
      }
    }
    else {
      profile[i] = raw_profile[i];
    }
  }
  mxFree(raw_profile);
  mxFree(multipliers);

}