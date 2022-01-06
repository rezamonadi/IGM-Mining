clear
fprintf('Setting paramters ...\n')
set_parameters
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
load(sprintf('%s/preloaded_qsos-mnp-%d', processed_directory(release),...
min_num_pixels));
fprintf('preparing voigt.c ...\n')

cd minFunc_2012
addpath(genpath(pwd))
mexAll
cd ..
% mex voigt.c -lcerf -I/rhome/reza66/bin/include/ -L/rhome/reza66/bin/lib/ % for hpcc
% mex voigt0.c -lcerf -I/rhome/reza66/bin/include/ -L/rhome/reza66/bin/lib/ % for hpcc
fprintf('Learning model ...\n')
lrn=0; % handle for learning or use saved learned model
training_set_name = 'mnp-200';
if lrn==1
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
    if(masking_CIV_region==1)
        train_ind = (filter_flags==0) & ismember(all_QSO_ID_dr7, half_ID);
    end

    learn_qso_model_dr7
else
    variables_to_load   = {'release', 'train_ind', 'max_noise_variance', ...
                            'minFunc_options', 'rest_wavelengths', 'mu', ...
                            'initial_M', 'M',  'log_likelihood', ...
                            'minFunc_output',   'restSN', 'coefficients',...
                            'latent', 'test_ind', 'prior_ind'};

    load(sprintf('%s/learned_model-%s.mat' , processed_directory(release),...
                                     training_set_name), variables_to_load{:});
end

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
load(sprintf('%s/preloaded_qsos', processed_directory(release)));
process_qsos_dr7



