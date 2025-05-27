% set_parameters: sets various parameters for the CIV detection 
% pipeline
% Desined for using DR16 spectra  
%flags for changes
extrapolate_subdla = 0; %0 = off, 1 = on
add_proximity_zone = 0;
integrate          = 1;
optTag = [num2str(integrate), num2str(extrapolate_subdla), num2str(add_proximity_zone)];

% physical constants
mgii_2796_wavelength = 2796.3543;		 % MgII transition wavelength  Å
mgii_2803_wavelength =  2803.5315; 		 % MgII transition wavelength  Å
speed_of_light = 299792458;                   % speed of light                     m s⁻¹

% converts relative velocity in km s^-1 to redshift difference
kms_to_z = @(kms) (kms * 1000) / speed_of_light;

% utility functions for redshifting
emitted_wavelengths = ...
    @(observed_wavelengths, z) (observed_wavelengths / (1 + z));

observed_wavelengths = ...
    @(emitted_wavelengths,  z) ( emitted_wavelengths * (1 + z));

%-------dr7
release = 'dr7';
% download Cooksey's dr7 spectra from this page: 
% http://www.guavanator.uhh.hawaii.edu/~kcooksey/SDSS/CIV/index.html 
% go to table: "SDSS spectra of the sightlines surveyed for C IV."
file_loader = @(mjd, plate, fiber_id) ...
  (read_spec_dr7(sprintf('data/dr7/spectro/1d_26/%04i/1d/spSpec-%05i-%04i-%03i.fit',...
  plate, mjd, plate, fiber_id)));


% % file loading parameters
loading_min_lambda = 1310;          % range of rest wavelengths to load  Å
loading_max_lambda = 2850;                    
% The maximum allowed is set so that even if the peak is redshifted off the end, the
% quasar still has data in the range

% preprocessing parameters
%z_qso_cut      = 2.15;                   % filter out QSOs with z less than this threshold
z_qso_cut      = 0.36;                      % according to Seyffert z>0.36                      
min_num_pixels = 100;                         % minimum number of non-masked pixels

% normalization parameters
% I use 1216 is basically because I want integer in my saved filenames%
% normalization_min_lambda = 1216 - 40;              % range of rest wavelengths to use   Å
%normalization_min_lambda = 1420; 
%normalization_max_lambda = 1216 + 40;              %   for flux normalization
%normalization_max_lambda = 1475; 

% NEW NORMALIZTION PARAMETERS
% range of rest wavelengths to use   Å
                                          % Continuum Fitting Redshift
                                          % Range
normalization_min_lambda_1 = 4150;   % Å  z < .6
normalization_max_lambda_1 = 4250;   % Å  z < .6

normalization_min_lambda_2 = 3020;   % Å  .6 < z < 1.0
normalization_max_lambda_2 = 3100;   % Å  .6 < z < 1.0

normalization_min_lambda_3 = 2180;   % Å  1.0 < z < 2.5   was 2150-2250
normalization_max_lambda_3 = 2250;   % Å  1.0 < z < 2.5

normalization_min_lambda_4 = 1420;   % Å  2.5 < z < 4.7
normalization_max_lambda_4 = 1500;   % Å  2.5 < z < 4.7

% null model parameters
min_lambda         =  loading_min_lambda+1;                   % range of rest wavelengths to       Å
max_lambda         = loading_max_lambda-1;                    %   model
dlambda            = 0.5;                    % separation of wavelength grid      Å
k                  = 20;                      % rank of non-diagonal contribution
max_noise_variance = 0.5^2;                   % maximum pixel noise allowed during model training
h                  = 2;                     % masking par to remove CIV region 

% optimization parameters
minFunc_options =               ...           % optimization options for model fitting
    struct('MaxIter',     10000, ...
           'MaxFunEvals', 10000);

% MgII model parameters: parameter samples (for Quasi-Monte Carlo)
nAVG               = 20;                     % number of points added between two 
                                            % observed wavelengths to make the Voigt finer
num_MgII_samples         = 20000;                  % number of parameter samples
alpha                    = 0.90;                    % weight of KDE component in mixture
uniform_min_log_nMgII     = 12.0;                   % range of column density samples    [cm⁻²]
uniform_max_log_nMgII     = 15.5;                   % from uniform distribution
fit_min_log_nMgII         = uniform_min_log_nMgII;                   % range of column density samples    [cm⁻²]
fit_max_log_nMgII         = 15.5;                   % from fit to log PDF


min_sigma                = 10e5;                   % cm/s -> b/sqrt(2) -> min Doppler par from Cooksey
max_sigma                = 125e5;                   % cm/s -> b/sqrt(2) -> max Doppler par from Cooksey

vCut                     = 767; %500;                    % maximum cut velocity for MgII system 
% model prior parameters
prior_z_qso_increase = kms_to_z(30000);       % use QSOs with z < (z_QSO + x) for prior


% instrumental broadening parameters
width = 3;                                    % width of Gaussian broadening (# pixels)
pixel_spacing = 1e-4;                         % wavelength spacing of pixels in dex
% DLA model parameters: absorber range and model
num_lines = 2;                                % number of members of CIV series to use

max_z_cut = kms_to_z(vCut);                   % max z_DLA = z_QSO - max_z_cut
max_z_MgII = @(z_qso, max_z_cut) ...         % determines maximum z_DLA to search
     z_qso - max_z_cut*(1+z_qso);
min_z_cut = kms_to_z(vCut);                   % min z_DLA = z_Ly∞ + min_z_cut
min_z_MgII = @(wavelengths, z_qso) ...         % determines minimum z_DLA to search
    max(min(wavelengths) / mgii_2796_wavelength - 1,                          ...
        observed_wavelengths(1310, z_qso) / mgii_2796_wavelength - 1);
train_ratio =0.99;
sample_name = sprintf("N-%d-%d-Sigma-%d-%d-Num-%d",floor(fit_min_log_nMgII*100),floor(100*fit_max_log_nMgII), min_sigma,max_sigma, num_MgII_samples);
training_set_name = 'Seyfert';
testing_set_name = '1%-Masking-All-samp-20k'
                                                
max_MgII = 10;
dv_mask = 750;%350; % (km/s)

                 
% base directory for all data   

base_directory = 'data';
% utility functions for identifying various directories
distfiles_directory = @(release) ...
   sprintf('%s/%s/distfiles', base_directory, release);

spectra_directory   = @(release)...
   sprintf('%s/%s/spectra', base_directory, release);

processed_directory = @(release) ...
   sprintf('%s/%s/processed', base_directory, release);

MgII_catalog_directory = @(name) ...
   sprintf('%s/MgII_catalogs/%s/processed', base_directory, name);
