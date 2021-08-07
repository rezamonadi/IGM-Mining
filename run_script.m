clc

fprintf('Setting paramters ...\n')
set_parameters
fprintf('Building catalogs ...\n')
build_catalog
fprintf('Preloading QSOs ...\n')

% preload_qsos
load(sprintf('data/dr7/processed/preloaded_qsos-dl-%d.mat', fix(dlambda*100)));
fprintf('preparing voigt.c ...\n')

cd minFunc_2012
addpath(genpath(pwd))
mexAll
cd ..
mex voigt.c -lcerf
fprintf('preparing testing, training, and prior indeces ...\n')

% f = load(sprintf('%s/filter_flags', processed_directory(training_release)));
% filter_flags= f.filter_flags;

% half_ID = randsample(all_QSO_ID, int32(train_ratio*numel(all_QSO_ID)));
% test_ind = (~ismember(all_QSO_ID, half_ID)) & (filter_flags==0) & (all_RATING~=2);

% prior_ind = ((ismember(all_QSO_ID, c4_QSO_ID)) & (filter_flags==0) & ...
%     (ismember(all_QSO_ID, half_ID)));
% train_ind = (~ismember(all_QSO_ID, c4_QSO_ID) & (filter_flags==0) & ...
% ismember(all_QSO_ID, half_ID) );

fprintf('Learning model (%s) ...\n', training_set_name)

% learn_qso_model
variables_to_load = {'training_release', 'train_ind', 'max_noise_variance', ...
                    'minFunc_options', 'rest_wavelengths', 'mu', ...
                     'initial_M', 'M',  'log_likelihood', 'test_ind',...
                     'prior_ind' };

load(sprintf('%s/learned_model-tr-%s', processed_directory(training_release),...
                training_set_name), variables_to_load{:});

fprintf('Generating samples for integrating out parameters in the model...\n')
% WrSampL2

variables_to_load = {'offset_z_samples', 'sigma_samples',...
                     'log_nciv_samples', 'nciv_samples'};
load(sprintf('%s/civ_samples_WR', processed_directory(training_release)), ...
     variables_to_load{:});
c4_catalog_name = 'cooksey';
fprintf('Processing (%s) ...\n', training_set_name)
process_qsos
% % % % variables_to_load = {'training_release', 'training_set_name', ...
% % % %     'c4_catalog_name', 'prior_ind', 'release', ...
% % % %     'test_set_name', 'test_ind', 'prior_z_qso_increase', ...
% % % %     'max_z_cut', 'num_lines', 'min_z_c4s', 'max_z_c4s', ...
% % % %     'log_priors_no_c4', 'log_priors_c4', ...
% % % %     'log_likelihoods_no_c4', 'sample_log_likelihoods_c4', ...
% % % %     'log_likelihoods_c4', 'log_posteriors_no_c4', ...
% % % %     'log_posteriors_c4', 'model_posteriors', 'p_no_c4', ...
% % % %     'p_c4', 'map_z_c4', 'map_N_c4'};

% % % % filename = sprintf('%s/b-sample-processed_qsos_%s', ...
% % % %     processed_directory(release), ...
% % % %     test_set_name);

% % % % load(filename, variables_to_load{:});




