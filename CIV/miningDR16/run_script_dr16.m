clear
fprintf('Setting paramters ...\n')
num_quasars = 100;
cataloging = 0;
preloading = 0;
sampling   = 0;
plotting   = 1;
processing = 1;
merging    = 0;
EWer       =0;
pltP       =0;
CredInt    =0;
dv_mask    = 250;
HPCC = 0;
voigtPrep = 0;
maskType = 1;
priorType = 1;
ind_S=1;
saving=1;
cores = 6;
set_parameters_dr16
training_set_name
fprintf('Building catalogs ...\n')
if cataloging == 1
    build_catalog_dr16
end
variables_to_load= {'all_plate_dr16', 'all_mjd_dr16', 'all_fiber_dr16', ...
'all_QSO_ID_dr16', 'all_RA_dr16', 'all_DEC_dr16', 'all_zqso_dr16','all_BAL_PROB'};
load(sprintf('%s/catalog', processed_directory(releaseTest)), ...
    variables_to_load{:});



fprintf('preparing voigt.c ...\n')
if voigtPrep == 1 
    cd minFunc_2012
    addpath(genpath(pwd))
    mexAll
    cd ..
    
    mex voigt_iP.c -lcerf
      
end



catDR7 = load(sprintf('%s/catalog', processed_directory(releasePrior)));
filter_flagsDR7 = load(sprintf('%s/filter_flags', processed_directory(releasePrior)), ...
'filter_flags');
prior_ind = (filter_flagsDR7.filter_flags==0); 
all_z_civ_C13 = catDR7.all_z_civ1;
all_REW_1548_DR7 = catDR7.all_EW1;
all_REW_1550_DR7 = catDR7.all_EW2;


fprintf('Learning model ...\n')

variables_to_load = { 'max_noise_variance', ...
                   'minFunc_options', 'rest_wavelengths', 'mu', ...
                    'initial_M', 'M',  'log_likelihood', ...
                    };

load(sprintf('%s/learned_model-%s', processed_directory(releasePrior),...
                                     training_set_name), variables_to_load{:});

fprintf('Generating samples for integrating out parameters in the model...\n')
% variables_to_load = {'offset_z_samples', 'offset_sigma_samples',...
%                      'log_nciv_samples', 'nciv_samples'};
if sampling==1
    generate_c4_samples
else
    load(sprintf('%s/civ_samples_N-%d-%d-Sigma-%d-%d-Num-%d.mat', processed_directory(releaseTest),...
        floor(fit_min_log_nciv*100), floor(fit_max_log_nciv*100), floor(min_sigma/1e5), floor(max_sigma/1e5), ...
        num_C4_samples));
end
    
fprintf(sprintf('%d Samples are generated\n', num_C4_samples));
% load preprocessed QSOs
fprintf('Preloading QSOs ...\n')
if preloading == 1
    preload_qsos_dr16
end
load(sprintf('%s/preloaded_qsos_%s', processed_directory(releaseTest), testing_set_name));


fprintf('preparing testing and prior indeces ...\n')
test_ind = (filter_flags==0);
if processing==1
    fprintf('processing QSO: %d to %d\n\n', ind_S, ind_S +  num_quasars-1);
    parpool('local', cores);
    process_qsos_dr16_par
    
end

if merging==1
    mergeProcessesDr16
end

if EWer==1
    EW_dr16_voigt
end

if pltP  ==1
   pltPost;
end

if CredInt==1
    CI_sort
end