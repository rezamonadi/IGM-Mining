
fprintf('Setting parameters ...\n')
% Processing is chunked so the script can be launched interactively or by a
% batch scheduler.  Adjust these three values for a batch job; `num_quasars`
% is clamped to the number of eligible preloaded spectra below.
% ind_S = 1;                   % plotting run: filtered test-set range
% num_quasars = 5000;           % covers all nine plotted_qso_ind entries
cores = 32;
cataloging  = 0;
preloading  = 0;
sampling    = 0;
plotting    = 0;
processing  = 1;
merging     = 0;
EWer        = 0;
pltP        = 0;
CredInt     = 0;
dv_mask     = 250;
voigtPrep = 0;
Averaging =1; 
maskType = 1;
priorType = 1;
saving=1;
set_parameters_dr16
training_set_name
fprintf('Building catalogs ...\n')
if cataloging == 1
    build_catalog_dr16
end
variables_to_load= {'all_plate_dr16', 'all_mjd_dr16', 'all_fiber_dr16', ...
'all_QSO_ID_dr16', 'all_RA_dr16', 'all_DEC_dr16', 'all_zqso_dr16','all_BAL_PROB'};
load('data/dr16/processed/catalogDR16.mat', variables_to_load{:});


fprintf('preparing voigt.c ...\n')
if voigtPrep == 1 
    cd minFunc_2012
    addpath(genpath(pwd))
    mexAll
    cd ..
    
    mex voigt_iP_transition.c -lcerf -I/rhome/rmona003/bin/include/ -L/rhome/rmona003/bin/lib64
      
end



catDR7 = load(sprintf('%s/catalog', processed_directory(releasePrior)));
filter_flagsDR7 = load('data/dr7/processed/filter_flags');
prior_ind = (filter_flagsDR7.filter_flags==0); 
all_z_civ_C13 = catDR7.all_z_civ1;
% all_REW_1548_DR7 = catDR7.all_EW1;
% all_REW_1550_DR7 = catDR7.all_EW2;


fprintf('Learning model ...\n')

variables_to_load = { 'max_noise_variance', ...
                   'minFunc_options', 'rest_wavelengths', 'mu', ...
                    'initial_M', 'M',  'log_likelihood', ...
                    };

load('data/dr7/processed/learned_model-C13-full.mat', variables_to_load{:});

fprintf('Generating samples for integrating out parameters in the model...\n')
% variables_to_load = {'offset_z_samples', 'offset_sigma_samples',...
%                      'log_nciv_samples', 'nciv_samples'};
if sampling==1
    generate_c4_samples
else
    load(sprintf('data/dr16/processed/civ_samples_N-%d-%d-Sigma-%d-%d-Num-%d.mat',...
        floor(fit_min_log_nciv*100), floor(fit_max_log_nciv*100), floor(min_sigma/1e5), floor(max_sigma/1e5), ...
        num_C4_samples));
end
    
fprintf(sprintf('%d Samples are generated\n', num_C4_samples));
% load preprocessed QSOs
fprintf('Preloading QSOs ...\n')
if preloading == 1
    preload_qsos_dr16_fixed
end
load('data/dr16/processed/preloaded_qsos_dr16.mat');


fprintf('preparing testing and prior indices ...\n')
test_ind = (filter_flags==0);

% `process_qsos_dr16_par` indexes the filtered test set.  Validate and clamp
% the requested chunk here so its saved filename always records the actual
% [start, end] range and the final chunk cannot run past the catalog.
all_num_quasars = nnz(test_ind);
if ind_S < 1 || ind_S > all_num_quasars
    error('run_script_dr16:InvalidStartIndex', ...
        'ind_S (%d) must be between 1 and %d eligible QSOs.', ...
        ind_S, all_num_quasars);
end
if num_quasars < 1 || floor(num_quasars) ~= num_quasars
    error('run_script_dr16:InvalidChunkSize', ...
        'num_quasars must be a positive integer.');
end
if cores < 1 || floor(cores) ~= cores
    error('run_script_dr16:InvalidCoreCount', ...
        'cores must be a positive integer.');
end
num_quasars = min(num_quasars, all_num_quasars - ind_S + 1);
ind_E = ind_S + num_quasars - 1;

% plotted_qso_ind = sort([393, 3187, 3387, 3246, 230, 710, 4616, 492, 1074]);
if processing==1
    fprintf('processing QSO: %d to %d\n\n', ind_S, ind_E);
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', cores);
    elseif pool.NumWorkers ~= cores
        warning('run_script_dr16:ExistingPool', ...
            ['Reusing the existing pool with %d workers (requested %d).'], ...
            pool.NumWorkers, cores);
    end
    if Averaging==1
        process_qsos_dr16_par_fixed
    else
        process_qsos_dr16_par_nofine
    end
    
end

if merging==1
    mergeProcessedDR16

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
