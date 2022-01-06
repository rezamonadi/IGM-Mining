clear
fprintf('Setting paramters ...\n')
set_parameters_dr7
training_set_name = 'dr12-p_c4-60-dl-25'
train_release = 'dr12'           
train_ratio = 0.8;
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
load(sprintf('%s/preloaded_qsos', processed_directory(train_release)));
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
% train_ind = (~ismember(all_QSO_ID_dr7, c4_QSO_ID) & (filter_flags==0) & ...
% ismember(all_QSO_ID_dr7, half_ID) ) & ind_RATING;
% if(masking_CIV_region==1)
%     train_ind = (filter_flags==0) & ismember(all_QSO_ID, half_ID);
% end

% load('data/dr7/processed/learned_model-norm-1420-1470-SN>4.mat', 'prior_ind', 'test_ind');


fprintf('Learning model ...\n')



% load preprocessed QSOs
% train_release = 'dr7'           
% preloaded_qsos_cat_name= sprintf('%s/preloaded_qsos',processed_directory(release));
% % learn_qso_model_dr7
variables_to_load   = {'release', 'train_ind', 'max_noise_variance', ...
                     'minFunc_options', 'rest_wavelengths', 'mu', ...
                     'initial_M', 'M',  'log_likelihood', ...
                     'minFunc_output',   'restSN', 'coefficients',...
                     'latent'};
% load(sprintf('%s/learned_model-%s' , processed_directory(train_release),...
%                                      training_set_name), variables_to_load{:});
       
load(sprintf('%s/learned_model-%s' , processed_directory(train_release),...
                                     training_set_name), variables_to_load{:});
release = 'dr7'  

% fprintf('Generating samples for integrating out parameters in the model...\n')
% % generate_c4_samples
if(RejectionSampling==0)
    variables_to_load = {'uniform_min_log_nciv', 'uniform_max_log_nciv', ...
                        'fit_min_log_nciv', 'fit_max_log_nciv', 'alpha', ...
                        'extrapolate_min_log_nciv', 'offset_z_samples',...
                        'offset_sigma_samples', 'log_nciv_samples', 'nciv_samples'};
    load(sprintf('%s/civ_samples', processed_directory(release)), variables_to_load{:});
else
    % WrSampL2_w_proposal_q
    variables_to_load = {'offset_z_samples', 'sigma_samples',...
    'log_nciv_samples', 'nciv_samples'};
    load(sprintf('%s/civ_samples_WR', processed_directory(release)), ...
                                                        variables_to_load{:});
end
fprintf('Processing ...\n')



Rsl =0; 
SnglMdl=0;
% Spline parameters
dv_civ_EW       = 700;   % km/s
dv_continuum_fit = 2000;  % km/s

testing_set_name = sprintf('SnglMdl-%d-Rsl-%d-Rjc-%d-nAvg-%d-N-%d-%d-sigma-%d-%d',...
                           SnglMdl, Rsl, RejectionSampling, nAVG, floor(uniform_min_log_nciv*100),...
                           floor(100*uniform_max_log_nciv), floor(min_sigma/1e5),...
                            floor(max_sigma/1e5));
load(sprintf('%s/preloaded_qsos', processed_directory(train_release)));
process_qsosSpline_dr7



