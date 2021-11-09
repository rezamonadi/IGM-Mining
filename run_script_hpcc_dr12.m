clear
fprintf('Setting paramters ...\n')
set_parameters
training_set_name
fprintf('Building catalogs ...\n')

build_catalog_dr12
variables_to_load= {'all_plate_dr12', 'all_mjd_dr12', 'all_fiber_dr12', ...
'all_QSO_ID_dr12', 'all_RA_dr12', 'all_DEC_dr12', 'all_zqso_dr12',};
load(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_load{:});
fprintf('Preloading QSOs ...\n')

preload_qsos_dr12


% load(sprintf('%s/preloaded_qsos_no_sky_mask', processed_directory(release)));
% fprintf('preparing voigt.c ...\n')

% cd minFunc_2012
% addpath(genpath(pwd))
% mexAll
% cd ..
% % mex voigt.c -lcerf -I/rhome/reza66/bin/include/ -L/rhome/reza66/bin/lib/ % for hpcc
% % mex voigt_dr.c -lcerf -I/rhome/reza66/bin/include/ -L/rhome/reza66/bin/lib/ % for hpcc
% % mex voigt0.c -lcerf -I/rhome/reza66/bin/include/ -L/rhome/reza66/bin/lib/ % for hpcc
% % mex voigt.c -lcerf
% % mex voigt_dr.c -lcerf
% fprintf('preparing testing, training, and prior indeces ...\n')

% % half_ID = randsample(all_QSO_ID, int32(train_ratio*numel(all_QSO_ID)));
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
% % % train_ind = (filter_flags==0) & ismember(all_QSO_ID, half_ID);

% fprintf('Learning model ...\n')
% % load preprocessed QSOs
% preloaded_qsos_cat_name= sprintf('%s/preloaded_qsos',processed_directory(release));
% % learn_qso_model
% variables_to_load = {'release', 'train_ind', 'max_noise_variance', ...
%                    'minFunc_options', 'rest_wavelengths', 'mu', ...
%                     'initial_M', 'M',  'log_likelihood', 'test_ind',...
%                     'prior_ind' };

% load(sprintf('%s/learned_model-%s' , processed_directory(release),...
%                                      training_set_name), variables_to_load{:});

% fprintf('Generating samples for integrating out parameters in the model...\n')
% % % generate_c4_samples
% variables_to_load = {'uniform_min_log_nciv', 'uniform_max_log_nciv', ...
%                     'fit_min_log_nciv', 'fit_max_log_nciv', 'alpha', ...
%                     'extrapolate_min_log_nciv', 'offset_z_samples',...
%                     'offset_sigma_samples', 'log_nciv_samples', 'nciv_samples'};
% load(sprintf('%s/civ_samples', processed_directory(release)), variables_to_load{:});

% % % % % WrSampL2_w_proposal_q
% % % % % variables_to_load = {'offset_z_samples', 'sigma_samples',...
% % % % % 'log_nciv_samples', 'nciv_samples'};
% % % % load(sprintf('%s/civ_samples_WR', processed_directory(training_release)), ...
% % % %                                                         variables_to_load{:});

% % % c4_catalog_name = 'cooksey';
% % % fprintf('Processing ...\n')
% % % % parpool('local', 32)


% % % % load('/home/reza/gpC4/data/dr7/processed/civ_samples-cooksey-sample-tr-50.mat', variables_to_load{:});
% % % % load redshifts from catalog to process

% % % % % load preprocessed QSOs
% % % % variables_to_load = {'all_wavelengths', 'all_flux', 'all_noise_variance', ...
% % % %     'all_pixel_mask', 'all_sigma_pixel'};
% % % % load(sprintf('%s/preloaded_qsos', processed_directory(release)), ...
% % % %     variables_to_load{:});


% % % % testing_set_name = sprintf('sigma-%d-%d-N-%d-%d', min_sigma/1e5, max_sigma/1e5, uniform_min_log_nciv*10 , uniform_max_log_nciv*10)
% testing_set_name = 'AVG-0-c4rmh-no-no-telluric-rm';
% % % % 
% process_qsosAvergaing





