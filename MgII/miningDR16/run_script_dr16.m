fprintf('Setting paramters ...\n')
num_quasars = 50;
cataloging  = 0;
preloading  = 0;
sampling    = 0;
plotting    = 1;
processing  = 1;
merging     = 0;
EWer        = 0;
pltP        = 0;
CredInt     = 0;
voigtPrep = 0;
maskType = 1;
priorType = 1;
saving=1;
ind_S =1;
nAVG = 20;
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
    
    mex voigt_iP.c -lcerf -I/rhome/rmona003/bin/include/ -L/rhome/rmona003/bin/lib64
      
end



catDR7 = load(sprintf('%s/catalog.mat', processed_directory(releasePrior)));
filter_flagsDR7 = load(sprintf('%s/filter_flags', processed_directory(releasePrior)), ...
'filter_flags');
prior_ind = (filter_flagsDR7.filter_flags==0); 

all_z_MgII_Seyfert = catDR7.all_z_MgII1;





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
    generate_MgII_samples
else
    load(sprintf('%s/MgII_samples_%s.mat', processed_directory(releaseTest), sample_name));
end
    
fprintf(sprintf('%d Samples are generated\n', num_MgII_samples));
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
    % parpool('local', cores);
    process_qsos_dr16_par
    
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
