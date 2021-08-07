clc
fprintf('Setting paramters ...\n')
set_parameters
train_ratio =0.5;
fprintf('Building catalogs ...\n')
build_catalog
fprintf('Preloading QSOs ...\n')

% preload_qsos
load('data/dr7/processed/preloaded_qsos.mat')
fprintf('preparing voigt.c ...\n')

cd minFunc_2012
addpath(genpath(pwd))
mexAll
cd ..
mex voigt.c -lcerf -I/rhome/reza66/bin/include/ -L/rhome/reza66/bin/lib/ % for hpcc
fprintf('preparing testing, training, and prior indeces ...\n')

f = load(sprintf('%s/filter_flags', processed_directory(training_release)));
filter_flags= f.filter_flags;

half_ID = randsample(all_QSO_ID, int32(train_ratio*numel(all_QSO_ID)));
test_ind = (~ismember(all_QSO_ID, half_ID)) & (filter_flags==0) & (all_RATING~=2);

prior_ind = ((ismember(all_QSO_ID, c4_QSO_ID)) & (filter_flags==0) & ...
    (ismember(all_QSO_ID, half_ID)));
train_ind = (~ismember(all_QSO_ID, c4_QSO_ID) & (filter_flags==0) & ...
ismember(all_QSO_ID, half_ID) );

fprintf('Learning model ...\n')

% learn_qso_model
variables_to_load = {'training_release', 'train_ind', 'max_noise_variance', ...
                   'minFunc_options', 'rest_wavelengths', 'mu', ...
                    'initial_M', 'M',  'log_likelihood', 'test_ind',...
                    'prior_ind' };

load(sprintf('%s/learned_model', processed_directory(training_release)),...
            variables_to_load{:});

fprintf('Generating samples for integrating out parameters in the model...\n')
generate_c4_samples


% variables_to_load = {'uniform_min_log_nciv', 'uniform_max_log_nciv', ...
%                     'fit_min_log_nciv', 'fit_max_log_nciv', 'alpha', ...
%                     'extrapolate_min_log_nciv', ...
%                     'offset_z_samples', 'offset_sigma_samples', 'log_nciv_samples', 'nciv_samples'};
% load(sprintf('%s/civ_samples-%s', processed_directory(training_release), training_set_name), ...
%     variables_to_load{:});
c4_catalog_name = 'cooksey';
fprintf('Processing ...\n')
parpool('local', 32)
process_qsos





