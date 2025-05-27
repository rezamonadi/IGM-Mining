clear
fprintf('Setting paramters ...\n')
num_C4_samples = 10000;
cataloging = 0;
preloading = 0;
learning   = 0;
sampling   = 0;
plotting   = 0;
processing = 0;
merging    = 0;
EWer       =0;
pltP       =0;
CredInt    =1;
dv_mask    = 300;
HPCC = 0;
voigtPrep = 1;
maskType = 1;
priorType = 1;
ind_S=1;
num_quasars =5000;
saving=0;
set_parameters_dr12
training_set_name
fprintf('Building catalogs ...\n')
if cataloging == 1
    build_catalog_dr12
end
variables_to_load= {'all_plate_dr12', 'all_mjd_dr12', 'all_fiber_dr12', ...
'all_QSO_ID_dr12', 'all_RA_dr12', 'all_DEC_dr12', 'all_zqso_dr12',};
load(sprintf('%s/catalog', processed_directory(releaseTest)), ...
    variables_to_load{:});



fprintf('preparing voigt.c ...\n')
if voigtPrep == 1 
    cd minFunc_2012
    addpath(genpath(pwd))
    mexAll
    cd ..
    
    mex voigt_iP.c -lcerf
    mex voigt_iP_1.c -lcerf
    mex voigt0.c -lcerf
    mex voigt1.c -lcerf
    
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

load(sprintf('%s/learned_model-%s', processed_directory(releaseTest),...
                                     training_set_name), variables_to_load{:});

fprintf('Generating samples for integrating out parameters in the model...\n')
% variables_to_load = {'offset_z_samples', 'offset_sigma_samples',...
%                      'log_nciv_samples', 'nciv_samples'};
load(sprintf('%s/civ_samples_N-%d-%d-Sigma-%d-%d-Num-%d.mat', processed_directory(releasePrior),...
    floor(fit_min_log_nciv*100), floor(fit_max_log_nciv*100), min_sigma, max_sigma, ...
    num_C4_samples));
    
fprintf(sprintf('%d Samples are generated\n', num_C4_samples));
% load preprocessed QSOs
fprintf('Preloading QSOs ...\n')
if preloading == 1
    preload_qsos_dr12
end
load(sprintf('%s/preloaded_qsos_%s', processed_directory(releaseTest), training_set_name));


fprintf('preparing testing and prior indeces ...\n')
test_ind = (filter_flags==0);
if processing==1
    fprintf('processing QSO: %d to %d\n\n', ind_S, ind_S +  num_quasars-1);
    % parpool('local', cores);
    process_qsos_dr12_par
    
end

if merging==1
    mergeProcessesDr12
end

if EWer==1
    EW_dr12_voigt
end

if pltP  ==1
   pltPost;
end

if CredInt==1
    CI_sort
end