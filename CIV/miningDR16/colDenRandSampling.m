%% C IV Absorber Sampling Script — Final Fixed Version (log-likelihoods)
close all; 
clear;
clc;
load('processed_qsos_dr16_parameters.mat');
% load('logN_grid.mat');  % Must be 1×10000 double
% fprintf('logN_grid size: %s\n', mat2str(size(logN_grid)));
load('civ_samples_N-1250-1610-Sigma-5-115-Num-10000.mat');

pr = load('preloaded_qsos_DR16.mat');
all_noise = pr.all_noise_variance;
all_noise = all_noise(savingCat.test_ind);
all_flux = pr.all_flux;
all_flux = all_flux(savingCat.test_ind);
SNR = zeros([nnz(savingCat.test_ind),1]);

for i=1:nnz(savingCat.test_ind)
    SNR(i) = median(all_flux{i}./sqrt(all_noise{i}));
end
m = matfile('processed_qsos_dr16_Likelihood.mat');
info = h5info(m.Properties.Source, '/LikelihoodCat/all_sample_log_likelihoods_c4L2');
[nqso_h5, nweights, max_abs_per_qso] = deal(info.Dataspace.Size(1), info.Dataspace.Size(2), info.Dataspace.Size(3));

logN_matrix = savingCat.all_map_N_c4L2;
z_abs_matrix = savingCat.all_map_z_c4L2;
p_matrix    = savingCat.all_p_c4;

sampledAbsorbers(size(logN_matrix,1)*7) = struct('qso_idx', [], 'SNR',[], 'abs_idx', [], 'logN_samples', [], 'p_c4', [], 'z_abs_samples', []);

fprintf('Sampling logN values for each absorber...\n');
count = 0;
nSamples =10;
%%
for qso_idx = 1:size(logN_matrix, 1)
    logN_row = logN_matrix(qso_idx, :);
    p_row    = p_matrix(qso_idx, :);
    z_abs_row = z_abs_matrix(qso_idx,:);
    valid_absorbers = find(~isnan(logN_row) & logN_row > 0 & z_abs_row>0);

    for a = valid_absorbers
        log_likelihoods = h5read(m.Properties.Source, ...
                '/LikelihoodCat/all_sample_log_likelihoods_c4L2', ...
                [qso_idx, 1, a], [1, nweights, 1]);
       
        log_likelihoods = double(log_likelihoods(:));

        
        % Convert log-likelihoods to safe probability weights
        weights = exp(log_likelihoods - max(log_likelihoods));
        weights_sum = sum(weights);
        
        weights = weights / weights_sum;

        
        idx_samples = randsample(1:nweights, nSamples, true, weights);
        

        logN_samples = log_nciv_samples(idx_samples);
        
      
        count = count + 1;
        sampledAbsorbers(count).qso_idx = qso_idx;
        sampledAbsorbers(count).abs_idx = a;
        sampledAbsorbers(count).logN_samples = logN_samples;
        sampledAbsorbers(count).SNR = SNR(qso_idx);
        sampledAbsorbers(count).p_c4 = p_row(a);
        sampledAbsorbers(count).z_abs_samples = z_abs_row(a);



      

        
    end

    if mod(qso_idx, 1000) == 0
        fprintf('Processed QSO %d/%d\n', qso_idx, size(logN_matrix, 1));
    end
end

%% Save results
save('sampled_logN_per_absorber_Smp10.mat', 'sampledAbsorbers', '-v7.3');
