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
savingCat.all_min_z_c4s                   = nan(all_num_quasars,1);
savingCat.all_max_z_c4s                   = nan(all_num_quasars,1);
savingCat.all_map_N_c4L2                  = nan(all_num_quasars, max_civ);
savingCat.all_map_sigma_c4L2              = nan(all_num_quasars, max_civ);
whos savingCat

filename = sprintf('%s/processedDR12/processed_qsos_tst_DR12_S_1_E_2000.mat', processed_directory(releaseTest));
loadCat = load(filename);
ind_S = 1
ind_E = 2000;   
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
savingCat.all_min_z_c4s(ind_S:ind_E,1)                      = loadCat.min_z_c4s;
savingCat.all_max_z_c4s(ind_S:ind_E,1)                      = loadCat.max_z_c4s;
savingCat.all_map_N_c4L2(ind_S:ind_E, :)                  = loadCat.map_N_c4L2;
savingCat.all_map_sigma_c4L2(ind_S:ind_E, :)              = loadCat.map_sigma_c4L2;
clear loadCat

filename = sprintf('%s/processedDR12/processed_qsos_tst_DR12_S_2001_E_5000.mat', processed_directory(releaseTest));
loadCat = load(filename);
ind_S = 2001
ind_E = 5000;        
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
savingCat.all_min_z_c4s(ind_S:ind_E,1)                      = loadCat.min_z_c4s;
savingCat.all_max_z_c4s(ind_S:ind_E,1)                      = loadCat.max_z_c4s;
savingCat.all_map_N_c4L2(ind_S:ind_E, :)                  = loadCat.map_N_c4L2;
savingCat.all_map_sigma_c4L2(ind_S:ind_E, :)              = loadCat.map_sigma_c4L2;
whos loadCat
clear loadCat

for ind_S = 5001:num_quasars:180001
    ind_S
    ind_E = num_quasars - 1 + ind_S;
    
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
    savingCat.all_min_z_c4s(ind_S:ind_E,1)                      = loadCat.min_z_c4s;
    savingCat.all_max_z_c4s(ind_S:ind_E,1)                      = loadCat.max_z_c4s;
    savingCat.all_map_N_c4L2(ind_S:ind_E, :)                  = loadCat.map_N_c4L2;
    savingCat.all_map_sigma_c4L2(ind_S:ind_E, :)              = loadCat.map_sigma_c4L2;
    clear loadCat
end

filename = sprintf('%s/processedDR12/processed_qsos_tst_DR12_S_185001_E_%d.mat',...
processed_directory(releaseTest), all_num_quasars);
loadCat = load(filename);
ind_S = 185001
ind_E = all_num_quasars;     
savingCat.all_p_c4(ind_S:ind_E, :)                        = loadCat.p_c4(1:ind_E - ind_S+1,:);
savingCat.all_p_c4L1(ind_S:ind_E, :)                      = loadCat.p_c4L1(1:ind_E - ind_S+1,:);
savingCat.all_p_no_c4(ind_S:ind_E, :)                     = loadCat.p_no_c4(1:ind_E - ind_S+1,:);
savingCat.all_log_posteriors_c4L2(ind_S:ind_E, :)         = loadCat.log_posteriors_c4L2(1:ind_E - ind_S+1,:);
savingCat.all_log_posteriors_c4L1(ind_S:ind_E, :)         = loadCat.log_posteriors_c4L1(1:ind_E - ind_S+1,:);
savingCat.all_log_posteriors_no_c4(ind_S:ind_E, :)        = loadCat.log_posteriors_no_c4(1:ind_E - ind_S+1,:);
savingCat.all_log_likelihoods_c4L2(ind_S:ind_E, :)        = loadCat.log_likelihoods_c4L2(1:ind_E - ind_S+1,:);
savingCat.all_sample_log_likelihoods_c4L2(ind_S:ind_E, :, :) = loadCat.sample_log_likelihoods_c4L2(1:ind_E - ind_S+1,:,:);
savingCat.all_log_likelihoods_no_c4(ind_S:ind_E, :)       = loadCat.log_likelihoods_no_c4(1:ind_E - ind_S+1,:);
savingCat.all_log_priors_c4(ind_S:ind_E, :)               = loadCat.log_priors_c4(1:ind_E - ind_S+1,:);
savingCat.all_log_priors_no_c4(ind_S:ind_E, :)            = loadCat.log_priors_no_c4(1:ind_E - ind_S+1,:);
savingCat.all_map_z_c4L2(ind_S:ind_E, :)                  = loadCat.map_z_c4L2(1:ind_E - ind_S+1,:);
savingCat.all_min_z_c4s(ind_S:ind_E,1)                      = loadCat.min_z_c4s(1:ind_E - ind_S+1);
savingCat.all_max_z_c4s(ind_S:ind_E,1)                      = loadCat.max_z_c4s(1:ind_E - ind_S+1);
savingCat.all_map_N_c4L2(ind_S:ind_E, :)                  = loadCat.map_N_c4L2(1:ind_E - ind_S+1,:);
savingCat.all_map_sigma_c4L2(ind_S:ind_E, :)              = loadCat.map_sigma_c4L2(1:ind_E - ind_S+1,:);   
savingCat.test_ind = test_ind;
clear loadCat

fprintf('Saving ... ')
save('processed_dr12.mat', 'savingCat', '-v7.3')


% num_quasars = nnz(test_ind);
% % Saving Kathy's data                
% zQSO_test = all_zqso(test_ind);

% % comparing Z and N
% Z_C13 = all_z_civ(test_ind,:);
% Rate_test = all_RATING(test_ind, :);
% dv= 350;
% tr=0.85;



% Full_catalog = ...
%     load(sprintf('%s/catalog', processed_directory(release)));

% all_wavelengths    =    all_wavelengths(test_ind);
% all_flux           =           all_flux(test_ind);
% all_pixel_mask     =     all_pixel_mask(test_ind);
% all_sigma_pixel    =    all_sigma_pixel(test_ind);
% z_qsos             =           all_zqso(test_ind);
% num_quasars        =                numel(z_qsos);

% % % preprocess model interpolants
% % % griddedInterpolant does an interpolation and gives a function handle
% % % based on the grided data. If the data is 2xD, {x,y} are the the same size as
% % % row  and columns. M is like M=f(x,y) and is like a matrix with each element
% % % M(i,j) = f(x(i), y(j))
% mu_interpolator = ...
%     griddedInterpolant(rest_wavelengths,        mu,        'linear');
% M_interpolator = ...
%     griddedInterpolant({rest_wavelengths, 1:k}, M,         'linear');

% % initialize results with nan
% REW_1548                        = nan(num_quasars, max_civ);
    
% for quasar_ind = 1:num_quasars
% %  for quasar_ind = 1:20
%     tic;
%     z_qso = z_qsos(quasar_ind);
%     fprintf('processing quasar %i/%i (z_QSO = %0.4f) ...\n', ...
%                               quasar_ind, num_quasars, z_qso);
    
%     this_wavelengths    =    all_wavelengths{quasar_ind};
%     this_wavelengths    =              this_wavelengths';
%     this_flux           =           all_flux{quasar_ind}; 
%     this_flux           =                     this_flux';
%     this_pixel_mask     =     all_pixel_mask{quasar_ind};
%     this_pixel_mask     =               this_pixel_mask';
%     this_sigma_pixel    =     all_sigma_pixel{quasar_ind};
%     this_sigma_pixel    =               this_sigma_pixel';
%     % 
%     % convert to QSO rest frame
%     this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qso);
    
%     unmasked_ind = (this_rest_wavelengths >= min_lambda) & ...
%         (this_rest_wavelengths <= max_lambda);% & (this_sigma_pixel>0);
%     % keep complete copy of equally spaced wavelengths for absorption
%     % computation
%     this_unmasked_wavelengths = this_wavelengths(unmasked_ind);
    
%     % [mask_ind] remove flux pixels with pixel_mask; pixel_mask is defined
%     % in read_spec_dr7.m
%     ind = unmasked_ind & (~this_pixel_mask);
%     this_wavelengths      =      this_wavelengths(ind);
%     this_rest_wavelengths = this_rest_wavelengths(ind);
%     this_flux             =             this_flux(ind);
%     this_sigma_pixel      =      this_sigma_pixel(ind);
    
%     % interpolate model onto given wavelengths
%     this_mu = mu_interpolator( this_rest_wavelengths);
   
    
    
%     padded_wavelengths = ...
%         [logspace(log10(min(this_unmasked_wavelengths)) - width * pixel_spaciprocessed_qsos_tst_DR12_S_1_2000.mat
%         this_unmasked_wavelengths;...
%         logspace(log10(max(this_unmasked_wavelengths)) + pixel_spacing,...
%         log10(max(this_unmasked_wavelengths)) + width * pixel_spacing,...
%         width)'...
%         ];

%     ind = (~this_pixel_mask(unmasked_ind));
%     lenW_unmasked = length(this_unmasked_wavelengths);
%     ind_not_remove = true(size(this_flux));
%     this_z_1548 = (this_wavelengths / civ_1548_wavelength) - 1;
%     for num_c4=1:max_civ
%         num_lines=1;
%         absorptionL1_fine= voigt_iP(finer(padded_wavelengths, nAVG),...
%             map_z_c4L2(quasar_ind,  num_c4), 10^map_N_c4L2(quasar_ind, num_c4),...
%             num_lines, map_sigma_c4L2(quasar_ind, num_c4), finer(this_sigma_pixel, nAVG));
%         absorptionL1 = Averager(absorptionL1_fine, nAVG, lenW_unmasked);
%         absorptionL1 = absorptionL1(ind);
%         % Equivalent width calculation 
%         this_rest_wavelengths = this_rest_wavelengths(ind);
%         dv_Doppler = sqrt(2)*3*map_sigma_c4L2(quasar_ind, num_c4);
%         dv_mid = (civ_1550_wavelength - civ_1548_wavelength)*0.5/civ_1548_wavelength*speed_of_light*100; % cm/s
%         dv_min = min(dv_mid, dv_Doppler);
%         indIntegration = abs(this_z_1548-map_z_c4L2(quasar_ind, num_c4))< dv_min; 
%         fprintf('size(W):%d-%d,size(abs):%d-%d', size(this_rest_wavelengths(indIntegration)),...
%         size(absorptionL1(indIntegration)));
%         REW_1548(quasar_ind, num_c4) = trapz(this_rest_wavelengths,...
%         1-absorptionL1);
%         fprintf('QSO:%d, c4:%d, z_civ:%.2f, p(CIV):%.2f, REW_1548:%.2f\n',...
%         quasar_ind, num_c4, map_z_c4L2(quasar_ind, num_c4), ...
%         p_c4(quasar_ind, num_c4),  REW_1548(quasar_ind, num_c4));
%     end
%     fprintf(' took %0.3fs.\n', toc);
% end
% save('REW_1548_DR7.mat', 'REW_1548', '-v7.3');
% % load('EW_DR7.mat', 'EW');



   
%     fig = figure('visible', 'off');
%     clf();
%     histogram(reshape(N_CIV_DR12, all_num_quasars*max_civ,1), 50)
%     set(gca,'YScale','log')
%     set(get(gca, 'XLabel'), 'String', 'N_{CIV}');
%     exportgraphics(fig, 'histDR12_N.png', 'Resolution', '1800')

%     fig = figure('visible', 'off');
%     clf();
%     histogram(reshape(Z_CIV_DR12, all_num_quasars*max_civ,1), 50)
%     set(gca,'YScale','log')
%     set(get(gca, 'XLabel'), 'String', 'Z_{CIV}');
%     exportgraphics(fig, 'histDR12_Z.png', 'Resolution', '1800')
%     fig = figure('visible', 'off');
%     clf();
%     histogram(reshape(sigma_CIV_DR12, all_num_quasars*max_civ,1), 50)
%     set(gca,'YScale','log')
%     set(get(gca, 'XLabel'), 'String', '\sigma_{CIV}');
%     exportgraphics(fig, 'histDR12_sigma.png', 'Resolution', '1800')
   
   

