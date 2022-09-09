% cataloging=0;preloading=0;sampling=0;learning=0;processing=1;plotting=1;
fprintf('Setting paramters ...\n')
set_parameters_dr7;
fprintf('testing set name: %s\ntraining set name: %s\n', testing_set_name,...
        training_set_name);

fprintf('Building catalogs ...\n')
if cataloging==1
    build_catalog_dr7;
end
variables_to_load= {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7', ...
'all_QSO_ID', 'all_RA', 'all_DEC', 'all_zqso',...
'all_N_civ','all_z_civ', 'all_RATING', 'c4_QSO_ID' };
load(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_load{:});
fprintf('Preloading QSOs ...\n')
if preloading==1
    preload_qsos_dr7
end
load(sprintf('%s/preloaded_qsos_%s.mat', processed_directory(release), training_set_name));

fprintf('preparing voigt.c ...\n')

if voigtPrep == 1 
    cd minFunc_2012
    addpath(genpath(pwd))
    mexAll
    cd ..
    if HPCC == 1
        mex voigt_iP.c -lcerf -I/rhome/rmona003/bin/include/ -L/rhome/rmona003/bin/lib64/ 
        % mex voigt0.c -lcerf -I/rhome/rmona003/bin/include/ -L/rhome/rmona003/bin/lib64/ 
    else 
        mex voigt_iP.c -lcerf
        % mex voigt0.c -lcerf
    end
end
fprintf('preparing testing, training, and prior indeces ...\n')
if learning==1
    half_ID = randsample(all_QSO_ID, int32(train_ratio*numel(all_QSO_ID)));
    test_ind = (~ismember(all_QSO_ID, half_ID)) & (filter_flags==0);
    prior_ind = ismember(all_QSO_ID, half_ID) & (filter_flags==0);
        
    % ind_RATING = filter_flags>10; % initializing with a all false array 
    % for i=1:length(ind_RATING)
        % if ~any(all_RATING(i,:)>=0)
            % ind_RATING(i)=true;
        % end
    % end
    train_ind = (~ismember(all_QSO_ID, c4_QSO_ID) & (filter_flags==0) & ...
    ismember(all_QSO_ID, half_ID) );

    fprintf('Learning model ...\n')
% load preprocessed QSOs
    preloaded_qsos_cat_name= sprintf('%s/preloaded_qsos_%s.mat',... 
                               processed_directory(release), training_set_name);
    learn_qso_model
end

if learning == 2 % --> full learning 
    fprintf('Learning model on the full catalog ...\n')
    preloaded_qsos_cat_name= sprintf('%s/preloaded_qsos_%s.mat',... 
                               processed_directory(release), training_set_name);
    learn_qso_model
end

variables_to_load = {'release', 'train_ind', 'max_noise_variance', ...
                   'minFunc_options', 'rest_wavelengths', 'mu', ...
                    'initial_M', 'M',  'log_likelihood', 'test_ind',...
                    'prior_ind' };

load(sprintf('%s/learned_model-%s.mat', processed_directory(release), ...
            training_set_name), variables_to_load{:});

fprintf('Generating samples for integrating out parameters in the model...\n');

if sampling==1
    generate_c4_samples
end

variables_to_load = {'uniform_min_log_nciv', 'uniform_max_log_nciv', ...
                    'fit_min_log_nciv', 'fit_max_log_nciv', 'alpha', ...
                    'extrapolate_min_log_nciv', 'offset_z_samples',...
                    'offset_sigma_samples', 'log_nciv_samples', 'nciv_samples'};

load(sprintf('%s/civ_samples_N_%d_%d_Sigma_%d_%d_num_%d.mat', processed_directory(release),...
uniform_min_log_nciv*100, uniform_max_log_nciv*100, min_sigma, max_sigma, ...
num_C4_samples), variables_to_load{:});

fprintf(sprintf('%d Samples are generated\n', num_C4_samples));

% load preprocessed QSOs
variables_to_load = {'all_wavelengths', 'all_flux', 'all_noise_variance', ...
    'all_pixel_mask', 'all_sigma_pixel'};
load(sprintf('%s/preloaded_qsos_%s.mat', processed_directory(release), training_set_name), ...
    variables_to_load{:});


% % % 
if processing==1
    parpool('local', cores);

    %--------------RUN 0 ------------
    % process_qsos_MF_multi_c4_1_single 
    % % MF code --> min_z_sepration = dv_mask, 1S model, multi C4 models,
    % % with averaging of 20, voigt0, no Occam razon for num_C4_samples, 
    % % prior multi--> Roman
    % % Purity and completeness was very terrible for this set-up:
    % % mx:1
    % % PM = 829
    % % GP = 339
    % % PM1GP1 = 285
    % % Purity = 84.07, Completness = 34.38
    % % mx:2
    % % PM = 829
    % % GP = 713
    % % PM1GP1 = 409
    % % Purity = 57.36, Completness = 49.34
    % % mx:3
    % % PM = 829
    % % GP = 1390
    % % PM1GP1 = 484
    % % Purity = 34.82, Completness = 58.38
    % % mx:4
    % % PM = 829
    % % GP = 2067
    % % PM1GP1 = 517
    % % Purity = 25.01, Completness = 62.36
    % % mx:5
    % % PM = 829
    % % GP = 2744
    % % PM1GP1 = 540
    % % Purity = 19.68, Completness = 65.14
    % % mx:6
    % % PM = 829
    % % GP = 3421
    % % PM1GP1 = 555
    % % Purity = 16.22, Completness = 66.95

    % -------------- RUN 1 --------------
    SingleLineModel = 1;
    plotting = 0; 

    process_qsos_dr7
    % -----Run-2-----------------------
    % voigtType = 0 -->  Roman Broadeing 
    % maskType = 0 --> mask S and D each run 
    % priorType = 0 --> Roman prior with same amount each run 
    % OccamRazor = 0 --> no MF Occan Razor for num_C4_civ

    % -----Run-3-----------------------
    % voigtType = 0 -->  Roman Broadeing 
    % maskType =  1 --> masking both most probable(D, S) and stop if Null is stronger 
    % priorType = 0 --> Roman prior with same amount each run 
    % OccamRazor = 0 --> no MF Occan Razor for num_C4_civ



    
  
    % Does the purity/completeness test...
end

if statTesting == 1
    statTest
end