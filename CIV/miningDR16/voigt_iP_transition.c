/*
 * voigt_iP_transition.c
 *
 * Instrument-broadened C IV Voigt profile for MATLAB MEX.
 *
 * Usage:
 *   profile = voigt_iP_transition(wavelengths, z, N, transition_mode, ...
 *                                 sigma, pixel_sigma, fine_factor)
 *
 * transition_mode:
 *   1  C IV 1548.1949 only
 *   2  C IV 1550.7700 only
 *   3  both transitions (the C IV doublet)
 *
 * wavelengths and pixel_sigma must be on the same fine grid.  pixel_sigma
 * (WDISP) is expressed in native SDSS-pixel units; fine_factor is the number
 * of fine samples per native pixel.  The input must be padded by
 * 3 * fine_factor samples at each end.  The output has that padding removed.
 */

#include "mex.h"
#include <cerf.h>
#include <math.h>

#define LAMBDAS_ARG          prhs[0]
#define Z_ARG                prhs[1]
#define N_ARG                prhs[2]
#define TRANSITION_MODE_ARG  prhs[3]
#define SIGMA_ARG            prhs[4]
#define PIXEL_SIGMA_ARG      prhs[5]
#define FINE_FACTOR_ARG      prhs[6]
#define PROFILE_ARG          plhs[0]

#define NUM_TRANSITIONS 2

static const double c = 2.99792458e10; /* cm s^-1 */

static const double transition_wavelengths[NUM_TRANSITIONS] = {
    1.5481949462890625e-5,
    1.55077001953125e-5
};

static const double leading_constants[NUM_TRANSITIONS] = {
    7.802895118381213e-08,
    3.899701297867750e-08
};

static const double gammas[NUM_TRANSITIONS] = {
    3.255002952981575e02,
    3.243136695286643e02
};

/* Native-pixel half-width of the truncated Gaussian LSF. */
static const int width = 3;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *lambdas;
    double *pixel_sigma;
    double *profile;
    double *raw_profile;
    double z, N, sigma, total, velocity;
    double pixel_width, kernel_q, kernel_q_squared, kernel_step;
    double kernel_weight, normalization, convolution;
    int transition_mode, first_transition, last_transition;
    int fine_factor, convolution_width;
    int i, j, ii;
    mwSize input_num_points, num_points, center;

    if (nrhs != 7) {
        mexErrMsgIdAndTxt("voigt_iP_transition:inputCount",
            "Expected 7 inputs: wavelengths, z, N, transition_mode, sigma, pixel_sigma, fine_factor.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("voigt_iP_transition:outputCount",
            "Expected one output.");
    }
    if (!mxIsDouble(LAMBDAS_ARG) || mxIsComplex(LAMBDAS_ARG) ||
        !mxIsDouble(PIXEL_SIGMA_ARG) || mxIsComplex(PIXEL_SIGMA_ARG)) {
        mexErrMsgIdAndTxt("voigt_iP_transition:type",
            "wavelengths and pixel_sigma must be real double arrays.");
    }

    lambdas = mxGetPr(LAMBDAS_ARG);
    pixel_sigma = mxGetPr(PIXEL_SIGMA_ARG);
    z = mxGetScalar(Z_ARG);
    N = mxGetScalar(N_ARG);
    sigma = mxGetScalar(SIGMA_ARG);
    transition_mode = (int)mxGetScalar(TRANSITION_MODE_ARG);
    fine_factor = (int)mxGetScalar(FINE_FACTOR_ARG);

    if (mxGetScalar(TRANSITION_MODE_ARG) != (double)transition_mode) {
        mexErrMsgIdAndTxt("voigt_iP_transition:transitionMode",
            "transition_mode must be an integer: 1 (1548), 2 (1550), or 3 (doublet).");
    }

    if (transition_mode == 1) {
        first_transition = 0;
        last_transition = 0;
    } else if (transition_mode == 2) {
        first_transition = 1;
        last_transition = 1;
    } else if (transition_mode == 3) {
        first_transition = 0;
        last_transition = 1;
    } else {
        mexErrMsgIdAndTxt("voigt_iP_transition:transitionMode",
            "transition_mode must be 1 (1548), 2 (1550), or 3 (doublet).");
    }

    if (fine_factor < 1 || mxGetScalar(FINE_FACTOR_ARG) != (double)fine_factor) {
        mexErrMsgIdAndTxt("voigt_iP_transition:fineFactor",
            "fine_factor must be a positive integer.");
    }

    input_num_points = mxGetNumberOfElements(LAMBDAS_ARG);
    if (mxGetNumberOfElements(PIXEL_SIGMA_ARG) != input_num_points) {
        mexErrMsgIdAndTxt("voigt_iP_transition:size",
            "pixel_sigma must have one value per fine-grid wavelength.");
    }

    convolution_width = width * fine_factor;
    if (input_num_points <= 2 * (mwSize)convolution_width) {
        mexErrMsgIdAndTxt("voigt_iP_transition:grid",
            "Fine-grid wavelength array is too short for LSF padding.");
    }

    num_points = input_num_points - 2 * (mwSize)convolution_width;
    PROFILE_ARG = mxCreateDoubleMatrix(num_points, 1, mxREAL);
    profile = mxGetPr(PROFILE_ARG);

    raw_profile = mxMalloc(input_num_points * sizeof(*raw_profile));
    for (i = 0; i < (int)input_num_points; i++) {
        total = 0.0;
        for (j = first_transition; j <= last_transition; j++) {
            velocity = lambdas[i] * c /
                (transition_wavelengths[j] * (1.0 + z)) / 1e8 - c;
            total -= leading_constants[j] * voigt(velocity, sigma, gammas[j]);
        }
        raw_profile[i] = exp(N * total);
    }

    for (i = 0; i < (int)num_points; i++) {
        center = (mwSize)i + (mwSize)convolution_width;

        if (!isfinite(pixel_sigma[center]) || pixel_sigma[center] <= 0.0) {
            profile[i] = raw_profile[center];
            continue;
        }

        /*
         * The Gaussian LSF is symmetric.  Let
         * q = exp(-0.5 / (fine_factor * WDISP)^2); then the weight at
         * fine-grid offset k is q^(k^2).  Generate successive weights with
         * q^(k^2) / q^((k-1)^2) = q^(2k-1), avoiding one exp call per tap.
         */
        pixel_width = fine_factor * pixel_sigma[center];
        kernel_q = exp(-0.5 / (pixel_width * pixel_width));
        kernel_q_squared = kernel_q * kernel_q;
        kernel_step = kernel_q;
        kernel_weight = 1.0;
        normalization = 1.0;
        convolution = raw_profile[center];

        for (ii = 1; ii <= convolution_width; ii++) {
            kernel_weight *= kernel_step;
            normalization += 2.0 * kernel_weight;
            convolution += kernel_weight *
                (raw_profile[center - (mwSize)ii] + raw_profile[center + (mwSize)ii]);
            kernel_step *= kernel_q_squared;
        }

        profile[i] = convolution / normalization;
    }

    mxFree(raw_profile);
}
