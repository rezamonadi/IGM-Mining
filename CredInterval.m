fprintf('loading...\n');
% fname = sprintf('%s/processed_qsos_dr12_N-1250-1610-S-35-115-nc-10k.mat', processed_directory(releaseTest));
% load(fname)  
    
% % loading full 3D poster
% sample_log_likelihoods_c4L2    = savingCat.all_sample_log_likelihoods_c4L2;
% map_z_c4L2     = savingCat.all_map_z_c4L2;
% map_N_c4L2     = savingCat.all_map_N_c4L2;
% map_sigma_c4L2 = savingCat.all_map_sigma_c4L2;  
% p_c4           = savingCat.all_p_c4; 
% z_qsos             =     all_zqso_dr12(test_ind);
% test_ind       = savingCat.test_ind;

% sample_log_likelihoods_c4L2 = sample_log_likelihoods_c4L2(1:10, :, :);
% map_z_c4L2 = map_z_c4L2(1:10,:);
% map_N_c4L2 = map_N_c4L2(1:10,:);
% map_sigma_c4L2 = map_sigma_c4L2(1:10,:);
% p_c4 = p_c4(1:10,:);
% z_qsos = z_qsos(1:10);
% vs = {'sample_log_likelihoods_c4L2', 'map_z_c4L2', 'map_N_c4L2','map_sigma_c4L2', 'p_c4', 'z_qsos'}
% save('short10.mat', vs{:});
load('short10.mat')
level = 0.95;
% num_quasars = nnz(test_ind);
num_quasars =10;
CI_Z = nan(num_quasars, 7, 3);
CI_N = nan(num_quasars, 7, 3);
CI_Sigma = nan(num_quasars, 7, 3);
CI_W_r_1548 = nan(num_quasars, 7, 3);
CI_W_r_1550 = nan(num_quasars, 7, 3);
sample_sigma_c4 = min_sigma + (max_sigma-min_sigma)*offset_sigma_samples;

CI_Z = nan(num_quasars, 7, 3);
CI_N = nan(num_quasars, 7, 3);
CI_Sigma = nan(num_quasars, 7, 3);
sample_sigma_c4 = min_sigma + (max_sigma-min_sigma)*offset_sigma_samples;
d_N = (15-12.88)/10000;
d_sigma = (max_sigma - min_sigma)/10000;
d_Wr = 5/100000;
for this_quasar_ind=1:num_quasars
    z_qso = z_qsos(this_quasar_ind);
    min_z_c4s(this_quasar_ind) = min_z_c4(1310, z_qso);
    max_z_c4s(this_quasar_ind) = max_z_c4(z_qso, max_z_cut);
    % binZ = (max_z_c4s(this_quasar_ind) - min_z_c4s(this_quasar_ind))/500;
    sample_z_c4 = ...
        min_z_c4s(this_quasar_ind) +  ...
        (max_z_c4s(this_quasar_ind) - min_z_c4s(this_quasar_ind)) * offset_z_samples;
    
    for num_c4=1:7
        if p_c4(this_quasar_ind, num_c4)>=0.0
            ind_CI_max = abs(sample_z_c4-map_z_c4L2(this_quasar_ind, num_c4))<kms_to_z(dv_mask)*(1+z_qso);
            d_z = range(sample_z_c4(ind_CI_max))/10000;
            P = sample_log_likelihoods_c4L2(this_quasar_ind, :, num_c4);
            P = P - max(P);
            P = exp(P);
            P = P(ind_CI_max);
            theta = [sample_z_c4(ind_CI_max)', sample_sigma_c4(ind_CI_max)',...
                log_nciv_samples(ind_CI_max)', W_r_1548_samples(ind_CI_max),...
                W_r_1550_samples(ind_CI_max)]; % parameter space
            dTheta = [d_z, d_sigma, d_N, d_Wr, d_Wr]; % search step for parameters
            [~,indLmax] = max(P);
            MAP_W_r_1548 = W_r_1548_samples(indLmax);
            MAP_W_r_1550 = W_r_1550_samples(indLmax);
            MAPtheta = [map_z_c4L2(this_quasar_ind, num_c4); ...
                        map_sigma_c4L2(this_quasar_ind, num_c4); ...
                        map_N_c4L2(this_quasar_ind, num_c4); MAP_W_r_1548; MAP_W_r_1550];
            CI = ciFunc(theta, dTheta, P, level, MAPtheta);

            CI_Z(this_quasar_ind,num_c4,1) = CI(1,1);
            CI_Z(this_quasar_ind,num_c4,2) = CI(1,2); 
            CI_Z(this_quasar_ind,num_c4,3) = CI(1,2) - CI(1,1);
            
            CI_Sigma(this_quasar_ind,num_c4,1) = CI(2,1);
            CI_Sigma(this_quasar_ind,num_c4,2) = CI(2,2);
            CI_Sigma(this_quasar_ind,num_c4,3) = CI(2,2) - CI(2,1);
            
            CI_N(this_quasar_ind,num_c4,2) = CI(3,2);
            CI_N(this_quasar_ind,num_c4,1) = CI(3,1);
            CI_N(this_quasar_ind,num_c4,3) = CI(3,2) - CI(3,1);

            CI_W_r_1548(this_quasar_ind,num_c4,2) = CI(4,2);
            CI_W_r_1548(this_quasar_ind,num_c4,1) = CI(4,1);
            CI_W_r_1548(this_quasar_ind,num_c4,3) = CI(4,2) - CI(4,1);

            CI_W_r_1550(this_quasar_ind,num_c4,2) = CI(5,2);
            CI_W_r_1550(this_quasar_ind,num_c4,1) = CI(5,1);
            CI_W_r_1550(this_quasar_ind,num_c4,3) = CI(5,2) - CI(5,1);

%             fprintf('QSO-%d-nc-%d\n', this_quasar_ind, num_c4);
%             fprintf('CI(z)=%f-%f\n', CI(1,:))
%             fprintf('CI(sigma)=%f-%f\n', CI(2,:))
%             fprintf('CI(N)=%f-%f\n', CI(3,:))
%             fprintf('iV=%d\n', CI(4,1))


            
        end
    end
end



vs = {'CI_Z', 'CI_Sigma', 'CI_N', 'CI_W_r_1548', 'CI_W_r_1550'};
save(sprintf('CredIntervals/CIs%d.mat', floor(100*level)), vs{:});


