set_parameters_dr12; dv_mask= 350; 
load(sprintf('%s/catalog', processed_directory('dr12')), 'all_zqso_dr12');
fprintf('loading...\n');
fname = 'data/dr12/processed/processed_qsos_tst_DR12_S_1_E_5000.mat';
% load(fname)  
variables_to_load = {'offset_z_samples', 'offset_sigma_samples',...
                     'log_nciv_samples', 'nciv_samples'};
load(sprintf('%s/civ_samples_N_%d_%d_Sigma_%d_%d_num_%d_alpha_%d.mat', processed_directory(releasePrior),...
    fit_min_log_nciv*100, fit_max_log_nciv*100, min_sigma, max_sigma, ...
    num_C4_samples, alpha*100), variables_to_load{:});
    
% loading full 3D poster
% sample_log_likelihoods_c4L2    = savingCat.all_sample_log_likelihoods_c4L2;
% map_z_c4L2     = savingCat.all_map_z_c4L2;
% map_N_c4L2     = savingCat.all_map_N_c4L2;
% map_sigma_c4L2 = savingCat.all_map_sigma_c4L2;  
% p_c4           = savingCat.all_p_c4; 
z_qsos             =     all_zqso_dr12(test_ind);
% test_ind       = savingCat.test_ind;
num_quasars = nnz(test_ind);
CI_Z = nan(5000, 7, 3);
CI_N = nan(5000, 7, 3);
CI_Sigma = nan(5000, 7, 3);
sample_sigma_c4 = min_sigma + (max_sigma-min_sigma)*offset_sigma_samples;
d_N = (15-12.88)/10000;
d_sigma = (max_sigma - min_sigma)/10000;
for this_quasar_ind=1:5000
    z_qso = z_qsos(this_quasar_ind);
    min_z_c4s(this_quasar_ind) = min_z_c4(1310, z_qso);
    max_z_c4s(this_quasar_ind) = max_z_c4(z_qso, max_z_cut);
    % binZ = (max_z_c4s(this_quasar_ind) - min_z_c4s(this_quasar_ind))/500;
    sample_z_c4 = ...
        min_z_c4s(this_quasar_ind) +  ...
        (max_z_c4s(this_quasar_ind) - min_z_c4s(this_quasar_ind)) * offset_z_samples;
    
    for num_c4=1:1
        if p_c4(this_quasar_ind, num_c4)>0.85
            ind_CI_max = abs(sample_z_c4-map_z_c4L2(this_quasar_ind, num_c4))<kms_to_z(dv_mask);
            d_z = range(sample_z_c4(ind_CI_max))/10000;
            P = sample_log_likelihoods_c4L2(this_quasar_ind, :, num_c4);
            P = P - max(P);
            P = exp(P);
            P = P(ind_CI_max);
            theta = [sample_z_c4(ind_CI_max)', sample_sigma_c4(ind_CI_max)', log_nciv_samples(ind_CI_max)']; % parameter space
            dTheta = [d_z, d_sigma, d_N]; % search step for parameters
            MAPtheta = [map_z_c4L2(this_quasar_ind, num_c4); ...
                        map_sigma_c4L2(this_quasar_ind, num_c4); ...
                        map_N_c4L2(this_quasar_ind, num_c4)];
            CI = CI3(theta, dTheta, P, 0.68, MAPtheta);

            % CI_Z(this_quasar_ind,num_c4,1) = min(sample_z_c4(inds_68 & ind_CI_max));
            % CI_Z(this_quasar_ind,num_c4,2) = max(sample_z_c4(inds_68 & ind_CI_max)); 
            % CI_Z(this_quasar_ind,num_c4,3) = CI_Z(this_quasar_ind,num_c4,2) - CI_Z(this_quasar_ind,num_c4,1);
            CI_Z(this_quasar_ind,num_c4,1) = CI(1,1);
            CI_Z(this_quasar_ind,num_c4,2) = CI(1,2); 
            CI_Z(this_quasar_ind,num_c4,3) = CI(1,2) - CI(1,1);
            
            % CI_Sigma(this_quasar_ind,num_c4,1) = min(sample_sigma_c4(inds_68 & ind_CI_max));
            % CI_Sigma(this_quasar_ind,num_c4,2) = max(sample_sigma_c4(inds_68 & ind_CI_max));
            % CI_Sigma(this_quasar_ind,num_c4,3) = CI_Sigma(this_quasar_ind,num_c4,2) - CI_Sigma(this_quasar_ind,num_c4,1);
            CI_Sigma(this_quasar_ind,num_c4,1) = CI(2,1);
            CI_Sigma(this_quasar_ind,num_c4,2) = CI(2,2);
            CI_Sigma(this_quasar_ind,num_c4,3) = CI(2,2) - CI(2,1);
            
            % CI_N(this_quasar_ind,num_c4,1) = min(log_nciv_samples(inds_68 & ind_CI_max));
            % CI_N(this_quasar_ind,num_c4,2) = max(log_nciv_samples(inds_68 & ind_CI_max));
            % CI_N(this_quasar_ind,num_c4,3) = CI_N(this_quasar_ind,num_c4,2) - CI_N(this_quasar_ind,num_c4,1);
            CI_N(this_quasar_ind,num_c4,2) = CI(3,2);
            CI_N(this_quasar_ind,num_c4,1) = CI(3,1);
            CI_N(this_quasar_ind,num_c4,3) = CI(3,2) - CI(3,1);

            fprintf('QSO-%d-nc-%d\n', this_quasar_ind, num_c4);
            fprintf('CI(z)=%f-%f\n', CI(1,:))
            fprintf('CI(sigma)=%f-%f\n', CI(2,:))
            fprintf('CI(N)=%f-%f\n', CI(3,:))
            fprintf('iV=%d\n', CI(4,1))


            
        end
    end
end


fig = figure();
allCI_Z = reshape(CI_Z(:,:,3), 5000*7,1);
histogram(allCI_Z)
set(get(gca, 'XLabel'), 'String', 'CI(Z_{CIV})');
exportgraphics(fig, 'testAllCIzCI3.png', 'Resolution', 800)

fig = figure();
allCI_Sigma = reshape(CI_Sigma(:,:,3), 5000*7,1);
histogram(allCI_Sigma)
set(get(gca, 'XLabel'), 'String', 'CI(\sigma)');
exportgraphics(fig, 'testAllCIsigmaCI3.png', 'Resolution', 800)

fig = figure();
allCI_N = reshape(CI_N(:,:,3), 5000*7,1);
histogram(allCI_N)
set(get(gca, 'XLabel'), 'String', 'CI(N)');
exportgraphics(fig, 'testAllCIncivCI3.png', 'Resolution', 800)






