z_qsos          = all_zqso_dr16(test_ind);
all_num_quasars = numel(z_qsos);
num_quasars     = 2000;
all_num_quasars = numel(z_qsos);


savingCat.all_p_c4                        = nan(all_num_quasars, max_civ);
savingCat.all_p_c4L1                      = nan(all_num_quasars, max_civ);
savingCat.all_p_no_c4                     = nan(all_num_quasars, max_civ);
savingCat.all_log_posteriors_c4L2         = nan(all_num_quasars, max_civ);
savingCat.all_log_posteriors_c4L1         = nan(all_num_quasars, max_civ);
savingCat.all_log_posteriors_no_c4        = nan(all_num_quasars, max_civ);
savingCat.all_log_likelihoods_c4L2        = nan(all_num_quasars, max_civ);
LikelihoodCat.all_sample_log_likelihoods_c4L2 = nan(all_num_quasars, num_C4_samples, max_civ);
savingCat.all_log_likelihoods_no_c4       = nan(all_num_quasars, max_civ);
savingCat.all_log_priors_c4               = nan(all_num_quasars, max_civ);
savingCat.all_log_priors_no_c4            = nan(all_num_quasars, max_civ);
savingCat.all_map_z_c4L2                  = nan(all_num_quasars, max_civ);
savingCat.all_map_z_c4L1                  = nan(all_num_quasars, max_civ);
savingCat.all_min_z_c4s                   = nan(all_num_quasars,1);
savingCat.all_max_z_c4s                   = nan(all_num_quasars,1);
savingCat.all_map_N_c4L2                  = nan(all_num_quasars, max_civ);
savingCat.all_map_N_c4L1                  = nan(all_num_quasars, max_civ);
savingCat.all_map_sigma_c4L2              = nan(all_num_quasars, max_civ);
savingCat.all_map_sigma_c4L1              = nan(all_num_quasars, max_civ);
savingCat.all_REW_1548_dr12               = nan(all_num_quasars, max_civ);

for ind_S = 1:num_quasars:250001
    ind_S
    ind_E = num_quasars - 1 + ind_S;
    if ind_S==250001
        ind_E = all_num_quasars;
    end
    
    filename = sprintf('%s/processed_qsos_tst_%s_S_%d_E_%d.mat',   processed_directory(releaseTest), testing_set_name, ind_S, ind_E);

    loadCat = load(filename);
    
    if ind_S~= 250001
        
        savingCat.all_p_c4(ind_S:ind_E, :)                        = loadCat.p_c4;
        savingCat.all_p_c4L1(ind_S:ind_E, :)                      = loadCat.p_c4L1;
        savingCat.all_p_no_c4(ind_S:ind_E, :)                     = loadCat.p_no_c4;
        savingCat.all_log_posteriors_c4L2(ind_S:ind_E, :)         = loadCat.log_posteriors_c4L2;
        savingCat.all_log_posteriors_c4L1(ind_S:ind_E, :)         = loadCat.log_posteriors_c4L1;
        savingCat.all_log_posteriors_no_c4(ind_S:ind_E, :)        = loadCat.log_posteriors_no_c4;
        savingCat.all_log_likelihoods_c4L2(ind_S:ind_E, :)        = loadCat.log_likelihoods_c4L2;
        LikelihoodCat.all_sample_log_likelihoods_c4L2(ind_S:ind_E, :, :) = loadCat.sample_log_likelihoods_c4L2;
        savingCat.all_log_likelihoods_no_c4(ind_S:ind_E, :)       = loadCat.log_likelihoods_no_c4;
        savingCat.all_log_priors_c4(ind_S:ind_E, :)               = loadCat.log_priors_c4;
        savingCat.all_log_priors_no_c4(ind_S:ind_E, :)            = loadCat.log_priors_no_c4;
        savingCat.all_map_z_c4L2(ind_S:ind_E, :)                  = loadCat.map_z_c4L2;
        savingCat.all_map_z_c4L1(ind_S:ind_E, :)                  = loadCat.map_z_c4L1;
        savingCat.all_min_z_c4s(ind_S:ind_E,1)                    = loadCat.min_z_c4s;
        savingCat.all_max_z_c4s(ind_S:ind_E,1)                    = loadCat.max_z_c4s;
        savingCat.all_map_N_c4L2(ind_S:ind_E, :)                  = loadCat.map_N_c4L2;
        savingCat.all_map_N_c4L1(ind_S:ind_E, :)                  = loadCat.map_N_c4L1;
        savingCat.all_map_sigma_c4L2(ind_S:ind_E, :)              = loadCat.map_sigma_c4L2;
        savingCat.all_map_sigma_c4L1(ind_S:ind_E, :)              = loadCat.map_sigma_c4L1;
        savingCat.test_ind                                        = loadCat.test_ind;
        % savingCat.all_REW_1548_dr12(ind_S:ind_E, :)               = loadCat.REW_1548_dr12;
    else
        num_quasars  = 983;

        savingCat.all_p_c4(ind_S:ind_E, :)                        = loadCat.p_c4(1:num_quasars, :);
        savingCat.all_p_c4L1(ind_S:ind_E, :)                      = loadCat.p_c4L1(1:num_quasars, :);
        savingCat.all_p_no_c4(ind_S:ind_E, :)                     = loadCat.p_no_c4(1:num_quasars, :);
        savingCat.all_log_posteriors_c4L2(ind_S:ind_E, :)         = loadCat.log_posteriors_c4L2(1:num_quasars, :);
        savingCat.all_log_posteriors_c4L1(ind_S:ind_E, :)         = loadCat.log_posteriors_c4L1(1:num_quasars, :);
        savingCat.all_log_posteriors_no_c4(ind_S:ind_E, :)        = loadCat.log_posteriors_no_c4(1:num_quasars, :);
        savingCat.all_log_likelihoods_c4L2(ind_S:ind_E, :)        = loadCat.log_likelihoods_c4L2(1:num_quasars, :);
        LikelihoodCat.all_sample_log_likelihoods_c4L2(ind_S:ind_E, :, :) = loadCat.sample_log_likelihoods_c4L2(1:num_quasars, :, :);
        savingCat.all_log_likelihoods_no_c4(ind_S:ind_E, :)       = loadCat.log_likelihoods_no_c4(1:num_quasars, :);
        savingCat.all_log_priors_c4(ind_S:ind_E, :)               = loadCat.log_priors_c4(1:num_quasars, :);
        savingCat.all_log_priors_no_c4(ind_S:ind_E, :)            = loadCat.log_priors_no_c4(1:num_quasars, :);
        savingCat.all_map_z_c4L2(ind_S:ind_E, :)                  = loadCat.map_z_c4L2(1:num_quasars, :);
        savingCat.all_map_z_c4L1(ind_S:ind_E, :)                  = loadCat.map_z_c4L1(1:num_quasars, :);
        savingCat.all_min_z_c4s(ind_S:ind_E,1)                    = loadCat.min_z_c4s(1:num_quasars,1);
        savingCat.all_max_z_c4s(ind_S:ind_E,1)                    = loadCat.max_z_c4s(1:num_quasars,1);
        savingCat.all_map_N_c4L2(ind_S:ind_E, :)                  = loadCat.map_N_c4L2(1:num_quasars, :);
        savingCat.all_map_N_c4L1(ind_S:ind_E, :)                  = loadCat.map_N_c4L1(1:num_quasars, :);
        savingCat.all_map_sigma_c4L2(ind_S:ind_E, :)              = loadCat.map_sigma_c4L2(1:num_quasars, :);
        savingCat.all_map_sigma_c4L1(ind_S:ind_E, :)              = loadCat.map_sigma_c4L1(1:num_quasars, :);
        savingCat.test_ind                                        = loadCat.test_ind;
    end
    clear loadCat
end



fprintf('Saving short Catalog ... \n')
fname = sprintf('%s/processed_qsos_dr16_parameters.mat', processed_directory(releaseTest));
save(fname, 'savingCat', '-v7.3')
clear savingCat

fprintf('Saving Likelihood Catalog ... \n')
fname = sprintf('%s/processed_qsos_dr16_Likelihood.mat', processed_directory(releaseTest));
save(fname, 'LikelihoodCat', '-v7.3')
clear LikelihoodCat

