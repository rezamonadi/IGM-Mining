clc
fprintf('Setting paramters ...\n')
set_parameters
train_ratio =0.8;
fprintf('Building catalogs ...\n')
build_catalog
fprintf('Preloading QSOs ...\n')

% preload_qsos
load(sprintf('%s/preloaded_qsos', processed_directory(release)));
fprintf('preparing voigt.c ...\n')

cd minFunc_2012
addpath(genpath(pwd))
mexAll
cd ..
mex voigt.c -lcerf -I/rhome/reza66/bin/include/ -L/rhome/reza66/bin/lib/ % for hpcc
% mex voigt.c -lcerf
fprintf('preparing testing, training, and prior indeces ...\n')

% half_ID = randsample(all_QSO_ID, int32(train_ratio*numel(all_QSO_ID)));
% test_ind = (~ismember(all_QSO_ID, half_ID)) & (filter_flags==0);

% prior_ind = ((ismember(all_QSO_ID, c4_QSO_ID)) & (filter_flags==0) & ...
%     (ismember(all_QSO_ID, half_ID)));
% train_ind = (~ismember(all_QSO_ID, c4_QSO_ID) & (filter_flags==0) & ...
% ismember(all_QSO_ID, half_ID) );

fprintf('Learning model ...\n')
% load preprocessed QSOs
preloaded_qsos_cat_name= sprintf('%s/preloaded_qsos',processed_directory(release));
% learn_qso_model
variables_to_load = {'training_release', 'train_ind', 'max_noise_variance', ...
                   'minFunc_options', 'rest_wavelengths', 'mu', ...
                    'initial_M', 'M',  'log_likelihood', 'test_ind',...
                    'prior_ind' };

load(sprintf('%s/learned_model-tr-%s', processed_directory(training_release), ...
        training_set_name), variables_to_load{:});

fprintf('Generating samples for integrating out parameters in the model...\n')
% generate_c4_samples
% variables_to_load = {'uniform_min_log_nciv', 'uniform_max_log_nciv', ...
%                     'fit_min_log_nciv', 'fit_max_log_nciv', 'alpha', ...
%                     'extrapolate_min_log_nciv', ...
%                     'offset_z_samples', 'offset_sigma_samples', 'log_nciv_samples', 'nciv_samples'};
% load(sprintf('%s/civ_samples-%s', processed_directory(training_release), training_set_name), ...
%     variables_to_load{:});

% WrSampL2_w_proposal_q
variables_to_load = {'offset_z_samples', 'sigma_samples',...
'log_nciv_samples', 'nciv_samples'};
load(sprintf('%s/civ_samples_WR', processed_directory(training_release)), ...
                                                        variables_to_load{:});

c4_catalog_name = 'cooksey';
fprintf('Processing ...\n')
% parpool('local', 32)


% load('/home/reza/gpC4/data/dr7/processed/civ_samples-cooksey-sample-tr-50.mat', variables_to_load{:});
% load redshifts from catalog to process
catalog = load(sprintf('%s/catalog', processed_directory(release)));

% load preprocessed QSOs
variables_to_load = {'all_wavelengths', 'all_flux', 'all_noise_variance', ...
    'all_pixel_mask', 'all_sigma_pixel'};
load(sprintf('%s/preloaded_qsos', processed_directory(release)), ...
    variables_to_load{:});

process_qsos





