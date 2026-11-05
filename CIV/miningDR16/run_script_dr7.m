fprintf('Setting parameters ...\n')
set_parameters_dr7_fixed;

num_C4_samples = 10000;
cataloging = 1;
preloading = 0;
learning   = 0;
sampling   = 0;
plotting   = 0;
processing = 0;
merging    = 0;
EWer       = 0;
pltP       = 0;
CredInt    = 0;
statTesting = 0;
dv_mask    = 350;
voigtPrep = 0;
maskType = 1;
priorType = 1;
ind_S=1;
saving=1;
fprintf('testing set name: %s\ntraining set name: %s\n', testing_set_name,...
      training_set_name);

fprintf('Building catalogs ...\n')
if cataloging==1
    build_catalog_dr7;
end
variables_to_load= {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7', ...
'all_QSO_ID', 'all_RA', 'all_DEC', 'all_zqso', 'all_EW1', 'all_errEW1', ...
'all_N_civ','all_z_civ1', 'all_z_civ2', 'all_z_civ3', 'all_RATING', 'c4_QSO_ID'};
load('data/dr7/processed/catalog.mat',   variables_to_load{:});
fprintf('Preloading QSOs ...\n')
if preloading==1
    preload_qsos_dr7_fixed
end
load('data/dr7/processed/preloaded_qsos_dr7.mat');

fprintf('preparing voigt.c ...\n')
cd minFunc_2012
addpath(genpath(pwd))
cd ..
if voigtPrep == 1 
    
    mexAll
    
    mex voigt_iP_fixed.c -lcerf -I/rhome/rmona003/bin/include/ -L/rhome/rmona003/bin/lib64/ 
end
fprintf('preparing testing, training, and prior indeces ...\n')
size(c4_QSO_ID)
size(all_QSO_ID)
if learning==1
    train_ind = (~ismember(all_QSO_ID(:), c4_QSO_ID(:)) & (filter_flags(:)==0));

    fprintf('Learning model ...\n')
    fprintf('Training QSOs: %d of %d\n', nnz(train_ind), numel(train_ind));

% load preprocessed QSOs
    preloaded_qsos_cat_name= sprintf('data/dr7/processed/preloaded_qsos_dr7.mat');
    learn_qso_model_fixed
end



% variables_to_load = {'release', 'train_ind', 'max_noise_variance', ...
%                    'minFunc_options', 'rest_wavelengths', 'mu', ...
%                     'initial_M', 'M',  'log_likelihood', 'test_ind',...
%                     'prior_ind' };

% load(sprintf('%s/learned_model-%s.mat', processed_directory(release), training_set_name), variables_to_load{:});

% fprintf('Generating samples for integrating out parameters in the model...\n');

% if sampling==1
%     generate_c4_samples
% end

% variables_to_load = {'offset_z_samples', 'offset_sigma_samples'
%                      'log_nciv_samples', 'nciv_samples'};

% load(sprintf('%s/civ_samples_%s.mat', processed_directory(release), sample_name), variables_to_load{:});

% fprintf(sprintf('%d Samples are generated\n', num_C4_samples));

% % load preprocessed QSOs
% variables_to_load = {'all_wavelengths', 'all_flux', 'all_noise_variance', ...
%     'all_pixel_mask', 'all_sigma_pixel'};
% load(sprintf('%s/preloaded_qsos_%s.mat', processed_directory(release), training_set_name), ...
%     variables_to_load{:});
% % 

% % % % 
% if processing==1
%     parpool('local', 6);
%     % SingleLineModel = 1;

%     process_qsos_dr7
% end

% if statTesting==1
% 	statTest
% end

% if EWer==1
%     EWer_dr7
% end
