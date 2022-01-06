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

test_ind_dr12 = (filter_flags==0);

%Training on low p_c4 spectra of dr12 
% load('data/dr12/processed/processed_qsos_RDR12.mat', 'p_c4');
% load('data/dr12/processed/all_pc4.mat', 'all_pc4');
% all_pc4 = all_pc4';
% mask_c4 = all_pc4<0.60;
% train_ind = (filter_flags==0);

%  Loading trained model from dr7
variables_to_load   = {'minFunc_options', 'rest_wavelengths', 'mu', ...
                            'initial_M', 'M',  'log_likelihood'};
training_set_name = 'mnp-200';
load(sprintf('%s/learned_model-%s.mat' , processed_directory('dr7'),...
                                     training_set_name), variables_to_load{:});
fprintf('Learning model\n')
% load preprocessed QSOs
% preloaded_qsos_cat_name= sprintf('%s/preloaded_qsos',processed_directory(release));
% z_qsos             =        all_zqso_dr12(train_ind);
% tic;
% learn_qso_model
% toc;
% variables_to_load = {'release', 'train_ind', 'max_noise_variance', ...
%                     'minFunc_options', 'rest_wavelengths', 'mu', ...
%                      'initial_M', 'M',  'log_likelihood', 'test_ind',...
%                      'prior_ind' };
%  training_set_name ='norm-1420-1470-SN>4';
%  load(sprintf('data/dr7/processed/learned_model-%s' , training_set_name),...
%   variables_to_load{:});

% fprintf('Generating samples for integrating out parameters in the model...\n')
% % % % % generate_c4_samples
variables_to_load = {'uniform_min_log_nciv', 'uniform_max_log_nciv', ...
                     'fit_min_log_nciv', 'fit_max_log_nciv', 'alpha', ...
                     'extrapolate_min_log_nciv', 'offset_z_samples',...
                     'offset_sigma_samples', 'log_nciv_samples', 'nciv_samples'};
load('data/dr7/processed/civ_samples', variables_to_load{:});


fprintf('Generating samples for integrating out parameters in the model...\n')
if(RejectionSampling==0)
    
    % generate_c4_samples
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


SnglMdl=1;
% eqWidth =0;
% Spline parameters
SplFit = 0;
dv_civ_EW       = 700;   % km/s
dv_continuum_fit = 2000;  % km/s
testing_set_name = sprintf('SnglMdl-%d-f_L1_sigma-%.2f-f_L1_N-%.2f-Rsl-%d-Rjc-%d-nAvg-%d-N-%.2fd-%.2f-sigma-%.2f-%.2f',...
                           SnglMdl, f_L1_sigma, f_L1_nciv, Rsl, RejectionSampling,...
                           nAVG, uniform_min_log_nciv, uniform_max_log_nciv,...
                           min_sigma/1e5, max_sigma/1e5)
load(sprintf('%s/preloaded_qsos', processed_directory('dr12')));
process_qsos_dr12







