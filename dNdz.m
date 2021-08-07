
% Comparing column density and redshift 
% from our MAP analysis to C13 measurements
clc
clear
set_parameters;
build_catalog;

variables_to_load = {'training_release', 'training_set_name', ...
    'c4_catalog_name', 'prior_ind', 'release', ...
    'test_ind', 'prior_z_qso_increase', ...
    'max_z_cut', 'min_z_c4s', 'max_z_c4s', ...
    'log_priors_no_c4', 'log_priors_c4', ...
    'log_likelihoods_no_c4', 'sample_log_likelihoods_c4L1', ...
    'sample_log_likelihoods_c4L2','log_likelihoods_c4L1', 'log_likelihoods_c4L2'...
    'log_posteriors_no_c4', 'log_posteriors_c4L1', 'log_posteriors_c4L1',...
    'model_posteriors', 'p_no_c4', ...
    'p_L1', 'map_z_c4L1', 'map_N_c4L1', 'map_sigma_c4L1', ...
    'map_z_c4L2', 'map_N_c4L2', 'map_sigma_c4L2'};

filename = sprintf('%s/processed_qsos_R%s.mat', ...
    processed_directory(release), ...
     training_set_name);
load(filename, variables_to_load{:});
% load QSO model from training release
variables_to_load = {'rest_wavelengths', 'mu', 'M'};
load(sprintf('%s/learned_model-tr-%s',...
processed_directory(training_release), training_set_name)...
,variables_to_load{:});

% load C4 samples from training release
variables_to_load = {'sigma_samples', 'offset_z_samples', 'log_nciv_samples', 'nciv_samples'};
load(sprintf('%s/civ_samples_WR', processed_directory(training_release)), ...
     variables_to_load{:});
     
% load preprocessed QSOs
variables_to_load = {'all_wavelengths', 'all_flux', 'all_noise_variance', ...
                     'all_pixel_mask', 'all_sigma_pixel'};
load(sprintf('%s/preloaded_qsos-dl-5', processed_directory(release)), ...
     variables_to_load{:});

% building samples-> Z and sigma
% sample_sigma_c4 = min_sigma + (max_sigma-min_sigma)*offset_sigma_samples;
sample_sigma_c4 = sigma_samples;

% enable processing specific QSOs via setting to_test_ind

test_ind = test_ind & filter_flags==0;
all_wavelengths    =    all_wavelengths(test_ind);
all_flux           =           all_flux(test_ind);
all_noise_variance = all_noise_variance(test_ind);
all_pixel_mask     =     all_pixel_mask(test_ind);


num_quasars = sum(test_ind);
ID = all_QSO_ID(test_ind);
Z_Cooksey_compare=[0];
NCIV_Cooksey_compare = [0];
jj=0;
for quasar_ind=1:num_quasars

    this_ID = ID{quasar_ind};
    this_systems = ismember(c4_QSO_ID, this_ID);
    % convert to QSO rest frame
    z_qsos = all_zqso(test_ind);
    if(sum(this_systems)>0)
        this_c4s = NCIV(this_systems);
        if (this_c4s>0)
            jj=jj+1;

            this_Zs  = Z_c4(this_systems);
            matched_Z_c4= this_Zs(abs(this_Zs -map_z_c4L2(quasar_ind))==min(abs(this_Zs-map_z_c4L2(quasar_ind))));
            matched_N_c4= this_c4s(abs(this_Zs -map_z_c4L2(quasar_ind))==min(abs(this_Zs-map_z_c4L2(quasar_ind))));
            NCIV_Cooksey_compare(jj) = matched_N_c4 - map_N_c4L2(quasar_ind);
            Z_Cooksey_compare(jj) = matched_Z_c4 - map_z_c4L2(quasar_ind);
        end
    end
end
fig=figure();
histogram(NCIV_Cooksey_compare, 50)
set(get(gca, 'XLabel'), 'String', 'N(C13) - N(MAP)');
exportgraphics(fig, sprintf('%s-N(C13)-N(MAP).pdf',...
                            training_set_name),...
                            'ContentType','vector')
fig=figure();
histogram(Z_Cooksey_compare, 50)
set(get(gca, 'XLabel'), 'String', 'z(C13) -z(MAP)');
exportgraphics(fig, sprintf('%s-z(C13)-z(MAP).pdf',...
                            training_set_name),'ContentType','vector')

