z_qsos             =     all_zqso_dr12(test_ind);
all_num_quasars    =                numel(z_qsos);

num_quasars = 5000;
all_num_quasars    =                numel(z_qsos);


savingCat.all_p_c4                        = nan(all_num_quasars, max_civ);
savingCat.all_p_c4L1                      = nan(all_num_quasars, max_civ);
savingCat.all_p_no_c4                     = nan(all_num_quasars, max_civ);
savingCat.all_log_posteriors_c4L2         = nan(all_num_quasars, max_civ);
savingCat.all_log_posteriors_c4L1         = nan(all_num_quasars, max_civ);
savingCat.all_log_posteriors_no_c4        = nan(all_num_quasars, max_civ);
savingCat.all_log_likelihoods_c4L2        = nan(all_num_quasars, max_civ);
savingCat.all_sample_log_likelihoods_c4L2 = nan(all_num_quasars, num_C4_samples, max_civ);
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
savingCat.all_REW_1548_dr12              = nan(all_num_quasars, max_civ);

for ind_S = 1:num_quasars:180001
    ind_S
    ind_E = num_quasars - 1 + ind_S;
    if ind_S==180001
        ind_E = 185425;
    end
    
    fname = sprintf('S_%d_E_%d', ind_S, ind_E);
    filename = sprintf('%s/processedDR12/processed_qsos_tst_DR12_%s.mat', processed_directory(releaseTest), fname);
    loadCat = load(filename);
    
    savingCat.all_p_c4(ind_S:ind_E, :)                        = loadCat.p_c4;
    savingCat.all_p_c4L1(ind_S:ind_E, :)                      = loadCat.p_c4L1;
    savingCat.all_p_no_c4(ind_S:ind_E, :)                     = loadCat.p_no_c4;
    savingCat.all_log_posteriors_c4L2(ind_S:ind_E, :)         = loadCat.log_posteriors_c4L2;
    savingCat.all_log_posteriors_c4L1(ind_S:ind_E, :)         = loadCat.log_posteriors_c4L1;
    savingCat.all_log_posteriors_no_c4(ind_S:ind_E, :)        = loadCat.log_posteriors_no_c4;
    savingCat.all_log_likelihoods_c4L2(ind_S:ind_E, :)        = loadCat.log_likelihoods_c4L2;
    savingCat.all_sample_log_likelihoods_c4L2(ind_S:ind_E, :, :) = loadCat.sample_log_likelihoods_c4L2;
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
    savingCat.all_REW_1548_dr12(ind_S:ind_E, :)               = loadCat.REW_1548_dr12;

    clear loadCat
end



fprintf('Saving ... ')
fname = sprintf('%s/processed_qsos_dr12_N-1250-1610-S-35-115-nc-10k.mat', processed_directory(releaseTest));
save(fname, 'savingCat', '-v7.3')
clear savingCat