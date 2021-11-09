
% Comparing column density and redshift 
% from our MAP analysis to C13 measurements
% plot_a_processed_qso.m : plot sample of spectrum with model dla
% and the sample likelihoods
%
% Usage:
% ----
% % first load the catalogues and spectra
% addpath dr16q/
% load_processed_catalogsset_parameters;
clc
clear
set_parameters;
build_catalog;
variables_to_load= {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7', ...
'all_QSO_ID', 'all_RA', 'all_DEC', 'all_zqso', 'EW1', 'EW2',...
'all_N_civ','all_z_civ', 'all_RATING', 'dla_QSO_ID','log_posteriors_dla',...
 'log_posteriors_no_dla' };
load(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_load{:});

% testing_set_name = sprintf('sigma-%d-%d-N-%d-%d', min_sigma/1e5, max_sigma/1e5, uniform_min_log_nciv*10 , uniform_max_log_nciv*10)
testing_set_name = 'sigma-15-65'

filename = sprintf('data/dr7/processed/processed_qsos_R%s.mat', testing_set_name);
% filename = '/home/reza/gpc/data/dr7/processed/processed_qsos_Rsigma-20-100-N-135-158-.mat'; 
% filename ='data/dr7/processed/processed_qsos_RRvoigt-10000Smp0-tr-80-20-vCutReza-5000.mat';
% filename ='data/dr7/processed/processed_qsos_Rprior-fixed.mat';
variables_to_load = {'training_release', 'training_set_name', ...
    'c4_catalog_name', 'prior_ind', 'release', ...
    'test_ind', 'prior_z_qso_increase', ...
    'max_z_cut', 'min_z_c4s', 'max_z_c4s', ...
    'log_priors_no_c4', 'log_priors_c4', ...
    'log_likelihoods_no_c4',  ...
    'sample_log_likelihoods_c4','log_likelihoods_c4'...
    'log_posteriors_no_c4', 'log_posteriors_c4', ...
    'model_posteriors', 'p_no_c4', 'p_c4', ...
    'map_z_c4', 'map_N_c4', 'map_sigma_c4', 'EqW1', 'EqW2', 'DoubletRatio'};

load(filename, variables_to_load{:});
load('data/C4_catalogs/Cooksey_C4_cat/processed/CIV-cat.mat','c4_QSO_ID','Z_c4','NCIV');

% % load QSO model from training release
% variables_to_load = {'rest_wavelengths', 'mu', 'M'};
% load(sprintf('%s/learned_model-tr-%s',...
% processed_directory(training_release), training_set_name)...
% ,variables_to_load{:});

% % load C4 samples from training release
% variables_to_load = {'sigma_samples', 'offset_z_samples', 'log_nciv_samples', 'nciv_samples'};
% load(sprintf('%s/civ_samples_WR', processed_directory(training_release)), ...
%      variables_to_load{:});
     

% variables_to_load = {'uniform_min_log_nciv', 'uniform_max_log_nciv', ...
%                     'fit_min_log_nciv', 'fit_max_log_nciv', 'alpha', ...
%                     'extrapolate_min_log_nciv', 'offset_z_samples',...
%                     'offset_sigma_samples', 'log_nciv_samples', 'nciv_samples'};
% load(sprintf('%s/civ_samples', processed_directory(training_release)), variables_to_load{:});
     
% % load preprocessed QSOs
% variables_to_load = {'all_wavelengths', 'all_flux', 'all_noise_variance', ...
%                      'all_pixel_mask', 'all_sigma_pixel'};
% load(sprintf('%s/preloaded_qsos-dl-5', processed_directory(release)), ...
%      variables_to_load{:});

% building samples-> Z and sigma
% sample_sigma_c4 = min_sigma + (max_sigma-min_sigma)*offset_sigma_samples;
% sample_sigma_c4 = sigma_samples;

% enable processing specific QSOs via setting to_test_ind




num_quasars = sum(test_ind);
ID = all_QSO_ID(test_ind);
Z_Cooksey_compare=[0];
NCIV_Cooksey_compare = [0];
jj=0;
for quasar_ind=1:num_quasars
    all_c4_NCIV_test=all_N_civ(test_ind,:);
    all_c4_Z_test=all_z_civ(test_ind, :);
    this_c4s = all_c4_NCIV_test(quasar_ind,:);
    this_Zs  = all_c4_Z_test(quasar_ind,:);
    matched_Z_c4= this_Zs(abs(this_Zs -map_z_c4(quasar_ind))==min(abs(this_Zs-map_z_c4(quasar_ind))));
    
    matched_N_c4= this_c4s(abs(this_Zs -map_z_c4(quasar_ind))==min(abs(this_Zs-map_z_c4(quasar_ind))));
    this_ID = ID{quasar_ind};
    num_this_systems = nnz(all_c4_NCIV_test(quasar_ind,:)>0);
    % convert to QSO rest frame
    
    if(num_this_systems>0)
        jj=jj+1;
        NCIV_Cooksey_compare(jj) = matched_N_c4(1) - map_N_c4(quasar_ind);
        Z_Cooksey_compare(jj) = matched_Z_c4(1) - map_z_c4(quasar_ind);
        
    end
end
fig=figure();
histogram(NCIV_Cooksey_compare, 15)
set(get(gca, 'XLabel'), 'String', 'N^{CIV}_{C13} - N^{CIV}_{MAP}');
set(get(gca, 'YLabel'), 'String', 'Frequency');
set(gca, 'FontSize',20)
exportgraphics(fig, sprintf('%s-N(C13)-N(MAP).pdf',...
                            training_set_name),...
                            'ContentType','vector')
fig=figure();
histogram(Z_Cooksey_compare, 15)
set(get(gca, 'XLabel'), 'String', 'z^{CVI}_{C13}-z^{CIV}_{MAP}');
set(get(gca, 'YLabel'), 'String', 'Frequency');
set(gca, 'FontSize',20)
exportgraphics(fig, sprintf('%s-z(C13)-z(MAP).pdf',...
                            training_set_name),'ContentType','vector')

