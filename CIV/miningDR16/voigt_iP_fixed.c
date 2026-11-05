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
#define FINE_FACTOR_ARG      prhs[6]
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
  double native_offset;
  int num_lines, i, j, ii, jj;
  int fine_factor, convolution_width;
  mwSize input_num_points, num_points, center;

  if (nrhs != 7) {
    mexErrMsgIdAndTxt("voigt_iP_fixed:inputCount",
                      "Expected 7 inputs, including fine_factor.");
  }
  if (nlhs > 1) {
    mexErrMsgIdAndTxt("voigt_iP_fixed:outputCount",
                      "Expected one output.");
  }
  if (!mxIsDouble(LAMBDAS_ARG) || mxIsComplex(LAMBDAS_ARG) ||
      !mxIsDouble(pixel_sigma_ARG) || mxIsComplex(pixel_sigma_ARG)) {
    mexErrMsgIdAndTxt("voigt_iP_fixed:type",
                      "Wavelengths and pixel_sigma must be real doubles.");
  }

  /* get input */
  lambdas     = mxGetPr(LAMBDAS_ARG);                /* wavelengths             Å             */
  z           = mxGetScalar(Z_ARG);                  /* redshift                dimensionless */
  N           = mxGetScalar(N_ARG);                  /* column density          cm⁻²          */
  sigma       = mxGetScalar(sigma_ARG);              /* Doppler broadening  paramberter   cm/s*/
  pixel_sigma = mxGetPr(pixel_sigma_ARG);             /* pixel sigma from spSpec                */

  fine_factor = (int)mxGetScalar(FINE_FACTOR_ARG);
  if (fine_factor < 1 || mxGetScalar(FINE_FACTOR_ARG) != (double)fine_factor) {
    mexErrMsgIdAndTxt("voigt_iP_fixed:fineFactor",
                      "fine_factor must be a positive integer.");
  }

  num_lines = (int)(mxGetScalar(NUM_LINES_ARG));
  input_num_points = mxGetNumberOfElements(LAMBDAS_ARG);
  if (mxGetNumberOfElements(pixel_sigma_ARG) != input_num_points) {
    mexErrMsgIdAndTxt("voigt_iP_fixed:size",
                      "pixel_sigma must have one value per fine-grid wavelength.");
  }
  if (num_lines < 1 || num_lines > NUM_LINES) {
    mexErrMsgIdAndTxt("voigt_iP_fixed:numLines",
                      "num_lines must be 1 or 2.");
  }

  convolution_width = width * fine_factor;
  if (input_num_points <= 2 * (mwSize)convolution_width) {
    mexErrMsgIdAndTxt("voigt_iP_fixed:grid",
                      "Fine-grid wavelength array is too short for LSF padding.");
  }

  /* initialize output */
  num_points = input_num_points - 2 * (mwSize)convolution_width;
  PROFILE_ARG = mxCreateDoubleMatrix(num_points, 1, mxREAL);
  profile = mxGetPr(PROFILE_ARG);                /* absorption profile      dimensionless */

  /* to hold the profile before instrumental broadening */
  raw_profile = mxMalloc(input_num_points * sizeof(double));
  multipliers = mxMalloc(num_lines * sizeof(double));
  instrument_profile = mxMalloc((2*convolution_width+1)* sizeof(double));

  for (i = 0; i < num_lines; i++)
    multipliers[i] = c / (transition_wavelengths[i] * (1 + z)) / 1e8;

  /* compute raw Voigt profile */
  for (i = 0; i < input_num_points; i++) {
    /* apply each absorption line */
    total = 0;
    for (j = 0; j < num_lines; j++) {
      /* velocity relative to transition wavelength */
      velocity = lambdas[i] * multipliers[j] - c;
      total += -leading_constants[j] * voigt(velocity, sigma, gammas[j]);
    }

    raw_profile[i] = exp(N * total);
  }

  /* instrumental broadening */
  for (i = 0; i < num_points; i++){
    center = (mwSize)i + (mwSize)convolution_width;
    if (!isfinite(pixel_sigma[center]) || pixel_sigma[center] <= 0) {
      profile[i] = raw_profile[center];
      continue;
    }

    total = 0;
    for (ii = -convolution_width, jj = 0;
         ii <= convolution_width; ii++, jj++) {
      /* ii is a fine-grid offset; WDISP is measured in native pixels. */
      native_offset = ii / (double)fine_factor;
      instrument_profile[jj] = exp(-0.5 * native_offset * native_offset /
                                   (pixel_sigma[center] * pixel_sigma[center]));
      total += instrument_profile[jj];
    }

    for (jj = 0; jj < 2 * convolution_width + 1; jj++)
      instrument_profile[jj] /= total;

    for (ii = -convolution_width, jj = 0;
         ii <= convolution_width; ii++, jj++) {
      profile[i] += raw_profile[(mwSize)((int)center + ii)] *
                    instrument_profile[jj];
    }
  }
  mxFree(raw_profile);
  mxFree(multipliers);
  mxFree(instrument_profile);

}
