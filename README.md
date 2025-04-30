# CIV Detection Pipeline for DR7 Spectra

This repository contains the implementation of a pipeline for detecting CIV absorption systems in quasar spectra from the Sloan Digital Sky Survey (SDSS) Data Release 7 (DR7). The pipeline is designed to preprocess, normalize, and analyze spectra to identify CIV systems using a combination of physical models and statistical techniques.

## Overview of `set_parameters_dr7.m`

The `set_parameters_dr7.m` file is a configuration script that defines key parameters, constants, and utility functions for the CIV detection pipeline. It is intended to be included in the main pipeline scripts and provides a centralized location for managing all configurable aspects of the analysis.

### Purpose

The script sets up:
- Physical constants and utility functions for redshift calculations.
- Parameters for file loading, preprocessing, and normalization.
- Model parameters for detecting CIV systems.
- Optimization settings for fitting models.
- Dataset configuration for training and testing.

### Key Sections

#### 1. **Flags for Changes**
This section defines flags that control optional features in the pipeline:
- `extrapolate_subdla`: Enables or disables extrapolation for sub-DLAs.
- `add_proximity_zone`: Toggles the inclusion of proximity zones.
- `integrate`: Determines whether integration is applied.

#### 2. **Physical Constants**
Defines fundamental constants used in the analysis:
- `civ_1548_wavelength` and `civ_1550_wavelength`: Wavelengths of the CIV doublet transitions.
- `speed_of_light`: Speed of light in meters per second.
- `kms_to_z`: A utility function to convert velocity (km/s) to redshift difference.

#### 3. **Utility Functions for Redshifting**
Provides helper functions for converting between emitted and observed wavelengths:
- `emitted_wavelengths`: Calculates emitted wavelengths given observed wavelengths and redshift.
- `observed_wavelengths`: Calculates observed wavelengths given emitted wavelengths and redshift.

#### 4. **Data Release Information**
Specifies the data release being used (`dr7`) and includes a file loader function for accessing spectra files.

#### 5. **File Loading Parameters**
Defines the wavelength range for loading spectra:
- `loading_min_lambda`: Minimum rest wavelength to load (Å).
- `loading_max_lambda`: Maximum rest wavelength to load (Å).

#### 6. **Preprocessing Parameters**
Sets thresholds for preprocessing:
- `z_qso_cut`: Minimum quasar redshift for inclusion in the analysis.
- `min_num_pixels`: Minimum number of non-masked pixels required for a spectrum.

#### 7. **Normalization Parameters**
Specifies the wavelength range for flux normalization:
- `normalization_min_lambda` and `normalization_max_lambda`: Define the rest wavelength range for normalization.

#### 8. **Null Model Parameters**
Defines parameters for modeling the null hypothesis (no CIV absorption):
- `min_lambda` and `max_lambda`: Wavelength range for modeling.
- `dlambda`: Wavelength grid spacing.
- `k`: Rank of the non-diagonal contribution in the model.
- `max_noise_variance`: Maximum allowable pixel noise variance.

#### 9. **Optimization Parameters**
Configures optimization settings for model fitting:
- `minFunc_options`: A structure defining the maximum number of iterations and function evaluations.

#### 10. **C4 Model Parameters**
Defines parameters for the CIV detection model:
- `nAVG`: Number of points added for finer Voigt profile.
- `num_C4_samples`: Number of parameter samples for the model.
- `alpha`: Weight of the KDE component in the mixture model.
- `uniform_min_log_nciv` and `uniform_max_log_nciv`: Column density range (log scale) for CIV systems.
- `min_sigma` and `max_sigma`: Doppler parameter range (cm/s).
- `vCut`: Maximum velocity cut for CIV systems (km/s).
- `prior_z_qso_increase`: Redshift increase for prior calculations.

#### 11. **Instrumental Broadening Parameters**
Defines parameters for simulating instrumental broadening:
- `width`: Width of Gaussian broadening in pixels.
- `pixel_spacing`: Wavelength spacing of pixels in logarithmic scale.

#### 12. **DLA Model Parameters**
Specifies parameters for modeling Damped Lyman-Alpha (DLA) systems:
- `num_lines`: Number of CIV series lines to use.
- `max_z_cut`, `max_z_c4`, and `min_z_c4`: Functions for calculating redshift cuts.

#### 13. **Dataset Parameters**
Configures the training and testing datasets:
- `train_ratio`: Proportion of data used for training.
- `sample_name`: Naming convention for sample datasets.
- `training_set_name` and `testing_set_name`: Names for training and testing datasets.

#### 14. **Base Directory and Utility Functions**
Defines directory paths for accessing data:
- `base_directory`: Root directory for data storage.
- `distfiles_directory`, `spectra_directory`, `processed_directory`, and `c4_catalog_directory`: Functions for constructing paths to specific data directories.

---

## Repository Structure

The repository is organized as follows:
