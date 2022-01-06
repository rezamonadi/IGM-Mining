fprintf('Setting paramters ...\n')
set_parameters_dr12
training_set_name
fprintf('Building catalogs ...\n')

% build_catalog_dr12
variables_to_load= {'all_plate_dr12', 'all_mjd_dr12', 'all_fiber_dr12', ...
'all_QSO_ID_dr12', 'all_RA_dr12', 'all_DEC_dr12', 'all_zqso_dr12',};
load(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_load{:});
fprintf('Preloading QSOs ...\n')

% preload_qsos_dr12


load(sprintf('%s/preloaded_qsos', processed_directory(release)));
% fprintf('preparing voigt.c ...\n')

cd minFunc_2012
addpath(genpath(pwd))
mexAll
cd ..
% % mex voigt.c -lcerf -I/rhome/reza66/bin/include/ -L/rhome/reza66/bin/lib/ % for hpcc
% % mex voigt_dr.c -lcerf -I/rhome/reza66/bin/include/ -L/rhome/reza66/bin/lib/ % for hpcc
% mex voigt0.c -lcerf -I/rhome/reza66/bin/include/ -L/rhome/reza66/bin/lib/ % for hpcc
% % mex voigt.c -lcerf
% % mex voigt_dr.c -lcerf
fprintf('preparing testing, training, and prior indeces ...\n')

test_ind_dr12 = filter_flags==0;
% % test_ind = (~ismember(all_QSO_ID, half_ID)) & (filter_flags==0);
% % prior_ind = ((ismember(all_QSO_ID, c4_QSO_ID)) & (filter_flags==0) & ...
% %     (ismember(all_QSO_ID, half_ID)));
% % ind_RATING = filter_flags>10; % initializing with a all false array 
% % for i=1:length(ind_RATING)
% %     if ~any(all_RATING(i,:)>=0)
% %         ind_RATING(i)=true;
% %     end
% % end
% % train_ind = (~ismember(all_QSO_ID, c4_QSO_ID) & (filter_flags==0) & ...
% % ismember(all_QSO_ID, half_ID) ) & ind_RATING;
% %%trainin inds by removing CIV absorbtion regions
% % train_ind = (filter_flags==0) & ismember(all_QSO_ID, half_ID);

%Training on low p_c4 spectra of dr12 
% load('data/dr12/processed/processed_qsos_RDR12.mat', 'p_c4');
% train_ind = (filter_flags==0) & (p_c4<0.5); 
fprintf('Learning model ...\n')
% load preprocessed QSOs
% preloaded_qsos_cat_name= sprintf('%s/preloaded_qsos',processed_directory(release));
% z_qsos             =        all_zqso_dr12(train_ind);
% learn_qso_model
variables_to_load = {'release', 'train_ind', 'max_noise_variance', ...
                    'minFunc_options', 'rest_wavelengths', 'mu', ...
                     'initial_M', 'M',  'log_likelihood', 'test_ind',...
                     'prior_ind' };
 training_set_name ='norm-1420-1470-SN>4';
 load(sprintf('data/dr7/processed/learned_model-%s' , training_set_name), variables_to_load{:});

fprintf('Generating samples for integrating out parameters in the model...\n')
% % % % generate_c4_samples
variables_to_load = {'uniform_min_log_nciv', 'uniform_max_log_nciv', ...
                     'fit_min_log_nciv', 'fit_max_log_nciv', 'alpha', ...
                     'extrapolate_min_log_nciv', 'offset_z_samples',...
                     'offset_sigma_samples', 'log_nciv_samples', 'nciv_samples'};
 load('data/dr7/processed/civ_samples', variables_to_load{:});

% % % % % % WrSampL2_w_proposal_q
% % % % % % variables_to_load = {'offset_z_samples', 'sigma_samples',...
% % % % % % 'log_nciv_samples', 'nciv_samples'};
% % % % % load(sprintf('%s/civ_samples_WR', processed_directory(training_release)), ...
% % % % %                                                         variables_to_load{:});

% % % % % % load preprocessed QSOs
% % % % % variables_to_load = {'all_wavelengths', 'all_flux', 'all_noise_variance', ...
% % % % %     'all_pixel_mask', 'all_sigma_pixel'};
% % % % % load(sprintf('%s/preloaded_qsos', processed_directory(release)), ...
% % % % %     variables_to_load{:});


% % % % % testing_set_name = sprintf('sigma-%d-%d-N-%d-%d', min_sigma/1e5, max_sigma/1e5, uniform_min_log_nciv*10 , uniform_max_log_nciv*10)
release = 'dr12';
num_quasars = nnz(test_ind_dr12);
ind_end = num_quasars;
testing_set_name = sprintf('DR12-flag-i-%d-f-%d', ind_start, ind_end);
fprintf('testing\nRelease:%s\ntesting set name: %s\n', release, testing_set_name);
% % % % 
% load C4 catalog
training_release = 'dr7';
Full_catalog = ...
     load(sprintf('%s/catalog', processed_directory(training_release)));

process_qsosAvergaing_dr12





