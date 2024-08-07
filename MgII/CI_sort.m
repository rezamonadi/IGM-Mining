fprintf('loading...\n');
fname = sprintf('%s/processed_qsos_dr12_N-1250-1610-S-35-115-nc-10k.mat', processed_directory(releaseTest));
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
% load('data/dr7/processed/civ_samples_N-1250-1610-Sigma-3500000-11500000-Num-10000.mat');
% load('REW_1548_DR12_fine.mat');
% load('short10.mat');
ConfidenceLevel = 0.68;
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


numColumnDensityBins = 9;
minEdgeColumnDensityBins = uniform_min_log_nciv ;
widthColumnDensityBins = (uniform_max_log_nciv - uniform_min_log_nciv)/numColumnDensityBins;
BinsColumnDensity = linspace(minEdgeColumnDensityBins,...
                            minEdgeColumnDensityBins + numColumnDensityBins*widthColumnDensityBins,...
                            numColumnDensityBins + 1);
midBinColumnDensity = BinsColumnDensity(2:end) - 0.5*widthColumnDensityBins;

i_hist = 0; 
AllCountColumnDensity = zeros([numColumnDensityBins, 1]);
for this_quasar_ind=1:num_quasars
    z_qso = z_qsos(this_quasar_ind);
    min_z_c4s(this_quasar_ind) = min_z_c4(1310, z_qso);
    max_z_c4s(this_quasar_ind) = max_z_c4(z_qso, max_z_cut);
    sample_z_c4 = ...
        min_z_c4s(this_quasar_ind) +  ...
        (max_z_c4s(this_quasar_ind) - min_z_c4s(this_quasar_ind)) * offset_z_samples;

    ThisCountColumnDensity = zeros([numColumnDensityBins, 1]);
    for num_c4=1:7
        if p_c4(this_quasar_ind, num_c4)>=0.0
            i_hist = i_hist + 1;
            indThisSystem = abs(sample_z_c4-map_z_c4L2(this_quasar_ind, num_c4))<kms_to_z(dv_mask)*(1+z_qso);

            max_log_likelihoodL2 = max(sample_log_likelihoods_c4L2(this_quasar_ind, indThisSystem, num_c4));
            sample_probabilitiesL2 = ...
                exp(sample_log_likelihoods_c4L2(this_quasar_ind, indThisSystem, num_c4)  ... 
                - max_log_likelihoodL2);

            Weights =  exp(sample_log_likelihoods_c4L2(this_quasar_ind, :, num_c4)  ... 
            - max_log_likelihoodL2);
            
            Weights(~indThisSystem) = 0;
            
            
            ThisCountColumnDensity(:) = SampleBinner(Weights,...
                                            log_nciv_samples,...
                                            minEdgeColumnDensityBins,...
                                            widthColumnDensityBins,...
                                            numColumnDensityBins,...
                                            10000);
            AllCountColumnDensity = AllCountColumnDensity + ThisCountColumnDensity;
            

            max_log_likelihoodL2 = max(sample_log_likelihoods_c4L2(this_quasar_ind, indThisSystem, num_c4));
            sample_probabilitiesL2 = ...
                exp(sample_log_likelihoods_c4L2(this_quasar_ind, indThisSystem, num_c4)  ... 
                - max_log_likelihoodL2);

            

            CI = ciFunc_sort(sample_probabilitiesL2,...
                             ConfidenceLevel,...
                             log_nciv_samples(indThisSystem),...
                             sample_sigma_c4(indThisSystem),...
                             sample_z_c4(indThisSystem),...
                             W_r_1548_samples(indThisSystem),...
                             W_r_1550_samples(indThisSystem));

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

            
        end
    end
end



vs = {'CI_Z', 'CI_Sigma', 'CI_N', 'CI_W_r_1548', 'CI_W_r_1550'};
save('CredIntervals/CIs65.mat', vs{:});


vs = {'midBinColumnDensity', 'AllCountColumnDensity'};

save('ColumnDensityDistData.mat',vs{:})