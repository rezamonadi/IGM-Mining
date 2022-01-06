clear
fprintf('Setting paramters ...\n')
set_parameters_dr7_lrn5
training_set_name
fprintf('Building catalogs ...\n')
% build_catalog_dr7
% ----------dr7----------------------------------------
variables_to_load = {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7', ...
'all_QSO_ID_dr7', 'all_RA_dr7', 'all_DEC_dr7', 'all_zqso_dr7',...
 'EW1', 'EW2', 'all_z_civ_c13', 'all_N_civ_c13', 'all_RATING_c13', 'c4_QSO_ID'};
load(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_load{:});

fprintf('Preloading QSOs ...\n')

% preload_qsos_dr7
load(sprintf('%s/preloaded_qsos', processed_directory(release)));
fprintf('preparing voigt.c ...\n')

cd minFunc_2012
addpath(genpath(pwd))
mexAll
cd ..
% mex voigt.c -lcerf -I/rhome/reza66/bin/include/ -L/rhome/reza66/bin/lib/ % for hpcc
% mex voigt_dr.c -lcerf -I/rhome/reza66/bin/include/ -L/rhome/reza66/bin/lib/ % for hpcc
% mex voigt0.c -lcerf -I/rhome/reza66/bin/include/ -L/rhome/reza66/bin/lib/ % for hpcc
% mex voigt.c -lcerf
% mex voigt_dr.c -lcerf
fprintf('preparing testing, training, and prior indeces ...\n')

half_ID = randsample(all_QSO_ID_dr7, int32(train_ratio*numel(all_QSO_ID_dr7)));
test_ind = (~ismember(all_QSO_ID_dr7, half_ID)) & (filter_flags==0);

prior_ind = ((ismember(all_QSO_ID_dr7, c4_QSO_ID)) & (filter_flags==0) & ...
    (ismember(all_QSO_ID_dr7, half_ID)));
ind_RATING = filter_flags>10; % initializing with a all false array 
for i=1:length(ind_RATING)
    if ~any(all_RATING_c13(i,:)>=0)
        ind_RATING(i)=true;
    end
end
train_ind = (~ismember(all_QSO_ID_dr7, c4_QSO_ID) & (filter_flags==0) & ...
ismember(all_QSO_ID_dr7, half_ID) ) & ind_RATING;
% train_ind = (filter_flags==0) & ismember(all_QSO_ID, half_ID);

% loading prior_ind from dr7 

% load('data/dr7/processed/learned_model-norm-1420-1470-SN>4.mat', 'prior_ind', 'test_ind');


train_release = 'dr7'           
% load(sprintf('%s/learned_model-%s' , processed_directory(train_release),...
                                    %  training_set_name), variables_to_load{:});

fprintf('Learning model ...\n')
% load preprocessed QSOs
% preloaded_qsos_cat_name= sprintf('%s/preloaded_qsos',processed_directory(release));
learn_qso_model_dr7
% variables_to_load = {'release', 'max_noise_variance', ...
%                    'minFunc_options', 'rest_wavelengths', 'mu', ...
%                     'initial_M', 'M',  'log_likelihood' };

% train_release = 'dr12'           
       
% load(sprintf('%s/learned_model-%s' , processed_directory(train_release),...
%                                      training_set_name), variables_to_load{:});
% release = 'dr7'  
% fprintf('Generating samples for integrating out parameters in the model...\n')
% % generate_c4_samples
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


% % % % load('/home/reza/gpC4/data/dr7/processed/civ_samples-cooksey-sample-tr-50.mat', variable  s_to_load{:});
% % % % load redshifts from catalog to process

% % % % testing_set_name = sprintf('sigma-%d-%d-N-%d-%d', min_sigma/1e5, max_sigma/1e5, uniform_min_log_nciv*10 , uniform_max_log_nciv*10)

% % % % 
% % dv_civ_EW       = 700;   % km/s
% % dv_continuum_fit = 2000;  % km/s
% % testing_set_name = sprintf('sp-1-1000-dvc-%d-dvci-%d', dv_continuum_fit, dv_civ_EW);
% release
% load(sprintf('%s/preloaded_qsos', processed_directory(release)));
% testing_set_name = 'C13-ratio-10'
% % process_qsosSpline_dr7
% process_qsosAvergaing_dr7





