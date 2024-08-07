% process_qsos: run CIV detection algorithm on specified objects

% Averaging     -> yes
% Single model  -> yes
% Multi CIV     -> yes
% Multi Singlet -> no
%
% removing CIV found region after each run 
% also removing the map Z_civ
% load C4 catalog
PM_catalog = ...                                   
    load(sprintf('%s/catalog', processed_directory(releasePrior)));
% train_ind -> those LOSs without any problem (filter_flag==0) that does not have
% Civ and useful for training null model (a model without Civ absorption line)
% prior_ind -> those LOSs with Civ absorption and in the half part of test
% test_ind -> second half without any filter flag for testing null and absorption model
% on the LOSs that we know have Civ or not. So, we asses our algorithm in this way.

if (ischar(prior_ind))
    prior_ind = eval(prior_ind);
end

% My prior_ind here is already those OK sight of lines that have CIV
prior.z_qsos  = PM_catalog.all_zqso(prior_ind);
prior.c4_ind = prior_ind;
prior.z_c4 = PM_catalog.all_z_civ(prior_ind);

% filter out CIVs from prior catalog corresponding to region of spectrum below
% Ly-alpha QSO rest. In Roman's code, for detecting DLAs, instead of
% Ly-alpha they have Lyman-limit
for i = size(prior.z_c4)
    if (observed_wavelengths(civ_1548_wavelength , prior.z_c4(i)) < ...
            observed_wavelengths(min_lambda, prior.z_qsos(i)))
        prior.c4_ind(i) = false;
    end
end
prior = rmfield(prior, 'z_c4');

% enable processing specific QSOs via setting to_test_ind
if (ischar(test_ind))
    test_ind = eval(test_ind);
end

all_wavelengths    =    all_wavelengths(test_ind);
all_flux           =           all_flux(test_ind);
all_noise_variance = all_noise_variance(test_ind);
all_pixel_mask     =     all_pixel_mask(test_ind);
all_sigma_pixel    =    all_sigma_pixel(test_ind);
z_qsos             =     all_zqso_dr12(test_ind);
all_num_quasars    =                numel(z_qsos);




% preprocess model interpolants
% griddedInterpolant does an interpolation and gives a function handle
% based on the grided data. If the data is 2xD, {x,y} are the the same size as
% row  and columns. M is like M=f(x,y) and is like a matrix with each element
% M(i,j) = f(x(i), y(j))
mu_interpolator = ...
    griddedInterpolant(rest_wavelengths,        mu,        'linear');
M_interpolator = ...
    griddedInterpolant({rest_wavelengths, 1:k}, M,         'linear');

% initialize results with nan
min_z_c4s                   = nan(num_quasars, 1);
max_z_c4s                   = nan(num_quasars, 1);
log_priors_no_c4            = nan(num_quasars, max_civ);
log_priors_c4               = nan(num_quasars, max_civ);
log_likelihoods_no_c4       = nan(num_quasars, max_civ);
sample_log_likelihoods_c4L1 = nan(num_quasars, num_C4_samples);
sample_log_likelihoods_c4L2 = nan(num_quasars, num_C4_samples, max_civ);
log_likelihoods_c4L1        = nan(num_quasars, max_civ);
log_likelihoods_c4L2        = nan(num_quasars, max_civ);
log_posteriors_no_c4        = nan(num_quasars, max_civ);
log_posteriors_c4L1         = nan(num_quasars, max_civ);
log_posteriors_c4L2         = nan(num_quasars, max_civ);
map_N_c4L1                  = nan(num_quasars, max_civ);
map_N_c4L2                  = nan(num_quasars, max_civ);
map_z_c4L1                  = nan(num_quasars, max_civ);
map_z_c4L2                  = nan(num_quasars, max_civ);
map_sigma_c4L1              = nan(num_quasars, max_civ);
map_sigma_c4L2              = nan(num_quasars, max_civ);
p_c4                        = nan(num_quasars, max_civ);
p_c4L1                      = nan(num_quasars, max_civ);
p_no_c4                     = nan(num_quasars, max_civ);
REW_1548                    = nan(num_quasars, max_civ);
num_pixel_civ               = nan(num_quasars, max_civ, 2);
sample_sigma_c4 = min_sigma + (max_sigma-min_sigma)*offset_sigma_samples;
% sample_sigma_c4 = sigma_samples;

% plt_count=0;
% FN_IDs = importdata('FN-list.csv');
% in_Kathy_FN_list = ismember(ID, FN_IDs);
% N_civ_test = all_N_civ(test_ind,:);
% z_PM_test = all_z_civ(test_ind,:);
z_PM_prior = all_z_civ_C13(prior_ind,:);
ind_E = all_num_quasars; ind_S=1;


this_quasar_ind = 0;
for all_quasar_ind = 1:all_num_quasars
    
    this_quasar_ind = this_quasar_ind + 1;
    tic;
    z_qso = z_qsos(all_quasar_ind);
    fprintf('processing quasar %i/%i (z_QSO = %0.4f) ...\n', ...
                              all_quasar_ind, num_quasars, z_qso);
    
    this_wavelengths    =    all_wavelengths{all_quasar_ind};
    % this_wavelengths    =              this_wavelengths';
    this_flux           =           all_flux{all_quasar_ind}; 
    % this_flux           =                     this_flux';
    this_noise_variance = all_noise_variance{all_quasar_ind};
    % this_noise_variance =           this_noise_variance';
    this_pixel_mask     =     all_pixel_mask{all_quasar_ind};
    % this_pixel_mask     =  this_pixel_mask';
    this_sigma_pixel   = all_sigma_pixel{all_quasar_ind};
    % this_sigma_pixel    =               this_sigma_pixel';


    % fprintf('b/f...');
    % fprintf('size(this_w)=%d-%d\n', size(this_wavelengths));
    % fprintf('size(this_flux)=%d-%d\n', size(this_flux));
    % fprintf('size(this_noise_variance)=%d-%d\n', size(this_noise_variance));
    % fprintf('size(this_sigma_pixel)=%d-%d\n', size(this_sigma_pixel));
    % fprintf('size(this_pixel_mask)=%d-%d\n', size(this_pixel_mask));
    % convert to QSO rest frame
    this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qso);
    
    unmasked_ind = (this_rest_wavelengths >= min_lambda) & ...
        (this_rest_wavelengths <= max_lambda) & (this_sigma_pixel>0);
    % keep complete copy of equally spaced wavelengths for absorption
    % computation
    this_unmasked_wavelengths = this_wavelengths(unmasked_ind);
    
    % [mask_ind] remove flux pixels with pixel_mask; pixel_mask is defined
    % in read_spec_dr7.m
    ind = unmasked_ind & (~this_pixel_mask);
    this_wavelengths      =      this_wavelengths(ind);
    this_rest_wavelengths = this_rest_wavelengths(ind);
    this_flux             =             this_flux(ind);
    this_noise_variance   =   this_noise_variance(ind);
    this_sigma_pixel      =      this_sigma_pixel(ind);
    

    % fprintf('after...');
    % fprintf('size(this_w)=%d-%d\n', size(this_wavelengths));
    % fprintf('size(this_flux)=%d-%d\n', size(this_flux));
    % fprintf('size(this_noise_variance)=%d-%d\n', size(this_noise_variance));
    % fprintf('size(this_sigma_pixel)=%d-%d\n', size(this_sigma_pixel));

    % c4 existence prior
    less_ind = (prior.z_qsos < (z_qso + prior_z_qso_increase));
    less_systems = z_PM_prior(less_ind,:);

    this_num_quasars = nnz(less_ind);
    this_p_c4(1) = nnz(less_systems(:,1)>0 )/this_num_quasars;  % P(at least 1 CIV)
    for i=2:max_civ
        this_p_c4(i) = nnz(less_systems(:,i)>0)/nnz(less_systems(:,i-1)>0); % P(at least n CIV | at least n-1 CIV)
        if (this_p_c4(i-1)==0)
        this_p_c4(i) = 0;
        end
    end

    fprintf('\n');
    for i = 1:max_civ
        if this_p_c4(i) > 0.90
            this_p_c4(i) = 0.90; % practically -Inf
        end
        log_priors_no_c4(this_quasar_ind, i) = log(1 - this_p_c4(i)); 
        fprintf(' ...     p(%i  CIVs | z_QSO)       : %0.3f\n', i, this_p_c4(i));
        fprintf(' ...     p(no CIV  | z_QSO)       : %0.3f\n', exp(log_priors_no_c4(this_quasar_ind, i)) );
    end
    log_priors_c4(this_quasar_ind,:) = log(this_p_c4(:));

    fprintf(' took %0.3fs.\n', toc);
   
end
% compute model posteriors in numerically safe manner



% save results
variables_to_save = { 'log_priors_no_c4', 'log_priors_c4'};

filename = sprintf('corrected Priors.mat');

save(filename, variables_to_save{:}, '-v7.3');