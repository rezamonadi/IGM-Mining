
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
    'map_z_c4L2', 'map_N_c4L2', 'map_sigma_c4L2', 'DoubletRatio', 'EqW1', 'EqW2'};

filename = sprintf('%s/processed_qsos_R%s', ...
    processed_directory(release), ...
    training_set_name);
processed=matfile(filename);
load('data/C4_catalogs/Cooksey_C4_cat/processed/CIV-cat.mat','c4_QSO_ID','Z_c4','NCIV');

% load QSO model from training release
variables_to_load = {'rest_wavelengths', 'mu', 'M'};
load(sprintf('%s/learned_model-tr-%s',...
processed_directory(processed.training_release), processed.training_set_name)...
,variables_to_load{:});

% load C4 samples from training release
variables_to_load = {'sigma_samples', 'offset_z_samples', 'log_nciv_samples', 'nciv_samples'};
load(sprintf('%s/civ_samples_WR', processed_directory(processed.training_release)), ...
     variables_to_load{:});
     
% load preprocessed QSOs
variables_to_load = {'all_wavelengths', 'all_flux', 'all_noise_variance', ...
                     'all_pixel_mask', 'all_sigma_pixel'};
load(sprintf('%s/preloaded_qsos', processed_directory(release)), ...
     variables_to_load{:});

% building samples-> Z and sigma
% sample_sigma_c4 = min_sigma + (max_sigma-min_sigma)*offset_sigma_samples;
sample_sigma_c4 = sigma_samples;




test_ind = processed.test_ind;
all_wavelengths    =    all_wavelengths(test_ind);
all_flux           =           all_flux(test_ind);
all_noise_variance = all_noise_variance(test_ind);
all_pixel_mask     =     all_pixel_mask(test_ind);
all_sigma_pixel    =    all_sigma_pixel(test_ind);

catalog = load(sprintf('%s/catalog', processed_directory(release)));
z_qsos = catalog.all_zqso(test_ind);
load(sprintf('%s/filter_flags', processed_directory(release)), ...
        'filter_flags');


p_L1 = processed.p_L1;
a = 1e-5;
p_c4s = 1- processed.p_no_c4 -a*p_L1;
sample_log_likelihoods_c4L2s = processed.sample_log_likelihoods_c4L2;
map_z_c4L2s = processed.map_z_c4L2;
map_N_c4L2s = processed.map_N_c4L2;
map_sigma_c4L2s = processed.map_sigma_c4L2;

map_z_c4L1s = processed.map_z_c4L1;
map_N_c4L1s = processed.map_N_c4L1;
map_sigma_c4L1s = processed.map_sigma_c4L1;

        % Testing Voigt profile 
dir = 'EW'; 
mkdir(dir);
count =0;
nqso = nnz(test_ind)

for quasar_ind=1:nqso
   
    % fprintf('quasar_ind:%d\n',quasar_ind);
    this_wavelengths    =    all_wavelengths{quasar_ind};
    this_wavelengths    =              this_wavelengths';
    this_flux           =           all_flux{quasar_ind}; 
    this_flux           =                     this_flux';
    this_noise_variance = all_noise_variance{quasar_ind};
    this_noise_variance =           this_noise_variance';
    this_pixel_mask     =     all_pixel_mask{quasar_ind};
    this_pixel_mask     =               this_pixel_mask';
    this_sigma_pixel    =     all_sigma_pixel{quasar_ind};
    this_sigma_pixel    =               this_sigma_pixel';
    % convert to QSO rest frame
  
    this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qsos(quasar_ind));
        
    unmasked_ind = (this_rest_wavelengths >= min_lambda) & ...
        (this_rest_wavelengths <= max_lambda) & (this_sigma_pixel>0);
    % keep complete copy of equally spaced wavelengths for absorption
    % computation
    this_unmasked_wavelengths = this_wavelengths(unmasked_ind);

    % [mask_ind] remove flux pixels with pixel_mask; pixel_mask is defined
    % in read_spec_dr7.m
    ind = unmasked_ind & (~this_pixel_mask) & (this_sigma_pixel>0);

    this_wavelengths      =      this_wavelengths(ind);
    this_rest_wavelengths = this_rest_wavelengths(ind);
    this_flux             =             this_flux(ind);
    this_noise_variance   =   this_noise_variance(ind);
    this_sigma_pixel      =      this_sigma_pixel(ind);

    % interpolate model onto given wavelengths
    mu_interpolator = ...
        griddedInterpolant(rest_wavelengths,        mu,        'linear');
    this_mu = mu_interpolator( this_rest_wavelengths);


    min_z_c4s = min_z_c4(this_wavelengths, z_qsos(quasar_ind));
    max_z_c4s = max_z_c4(z_qsos(quasar_ind), max_z_cut);
    sample_z_c4 = min_z_c4s + (max_z_c4s - min_z_c4s) * offset_z_samples;

    % ensure enough pixels are on either side for convolving with
    % instrument profile
    padded_wavelengths = ...
        [logspace(log10(min(this_unmasked_wavelengths)) - width * pixel_spacing, ...
        log10(min(this_unmasked_wavelengths)) - pixel_spacing,         ...
        width)';                                                       ...
        this_unmasked_wavelengths;                                              ...
        logspace(log10(max(this_unmasked_wavelengths)) + pixel_spacing,         ...
        log10(max(this_unmasked_wavelengths)) + width * pixel_spacing, ...
        width)'                                                        ...
        ];
    
    % [mask_ind] to retain only unmasked pixels from computed absorption profile
    % this has to be done by using the unmasked_ind which has not yet
    % been applied this_pixel_mask.
    ind = (~this_pixel_mask(unmasked_ind));
    ID = all_QSO_ID(test_ind);
    this_ID = ID{quasar_ind};
    this_systems = ismember(c4_QSO_ID, this_ID);
    
    % sample_log_likelihoods_c4 = sample_log_likelihoods_c4L2s(quasar_ind);
    num_systems = sum(this_systems);
    if(processed.EqW1(quasar_ind,1)<1 | processed.EqW2(quasar_ind,1)<1) % FP
        % if(num_systems==0 & p_c4s(quasar_ind)<0.2 & filter_flags(quasar_ind)==0) % TN
        count=count+1;
        fprintf('EW:%d\n',count);
        fig = figure('visible', 'off', 'Position', [0,0,2024,1800]);
        clf();
        % subplot('position', [0.05 0.49 0.90 0.45]);
% construct dla_mu_map
        num_lines=2;
        % absorptionL2 = voigt(padded_wavelengths, map_z_c4L2s(quasar_ind), ...
        %     10^map_N_c4L2s(quasar_ind), num_lines, map_sigma_c4L2s(quasar_ind), ...
        %     this_sigma_pixel);
        L1 = voigt0(padded_wavelengths, map_z_c4L2s(quasar_ind), ...
                    10^map_N_c4L2s(quasar_ind), 1, map_sigma_c4L2s(quasar_ind));
        L2 = voigt1(padded_wavelengths, map_z_c4L2s(quasar_ind), ...
        10^map_N_c4L2s(quasar_ind), 1, map_sigma_c4L2s(quasar_ind));
        absorptionL2 = L1.*L2;
        absorptionL2 = absorptionL2(ind);
        c4_muL2 = this_mu.* absorptionL2;
               
        this_z_c4 = (this_wavelengths / 1550) - 1;
        subplot(2,1,1);
        xlabel('(observed wavelengths $\lambda$ (\AA) / 1549 (\AA)) - 1', 'FontSize', 14, 'Interpreter','latex');
        ylabel('normalized flux $\mathbf{y}$',                            'FontSize', 14, 'Interpreter','latex');
        
        
        plot(this_z_c4, absorptionL2, 'Color', 'g');
        % hold on
        % plot(this_z_c4, this_flux, 'Color', 'r');
        % hold on
        % legend('L2', 'Raw')
       
        xlim([min(this_z_c4), max(this_z_c4)])

        xlabel('Z')
        ylabel('Normalized Flux')
        tit = sprintf('ID:%s\nzqso=%.2f, EW(1548):%.2f, EW(1550):%.2f, DR:%.2f, sigma:%.2fkm/s, N:%.2f, P:%.2f',...
                     ID{quasar_ind}, z_qsos(quasar_ind), processed.EqW1(quasar_ind,1),...
                     processed.EqW2(quasar_ind,1), processed.DoubletRatio(quasar_ind,1),...
                     processed.map_sigma_c4L2(quasar_ind,1)/1e5, processed.map_N_c4L2(quasar_ind,1),...
                     p_c4s(quasar_ind));
        title(tit)
        title(sprintf('EW1:%.2e, EW2:%.2e, intEW:%.2f',...
             processed.EqW1(quasar_ind,1), processed.EqW2(quasar_ind,1), ...
             trapz(this_unmasked_wavelengths, 1-absorptionL2)))
        hold on 



        subplot(2,1,2);
        hold on
        norm_sample_log_likelihoods = processed.sample_log_likelihoods_c4L2(quasar_ind, :);
        norm_sample_log_likelihoods = norm_sample_log_likelihoods - max(norm_sample_log_likelihoods);
        norm_sample_log_likelihoods = norm_sample_log_likelihoods - log(sum(exp(norm_sample_log_likelihoods)));

        % s=scatter(sample_sigma_c4, log_nciv_samples,20,...
                    %  norm_sample_log_likelihoods, 'filled', 'DisplayName', 'sample likelihoods');
        s=scatter(sample_z_c4, log_nciv_samples,20,...
                     norm_sample_log_likelihoods, 'filled', 'DisplayName', 'sample likelihoods');
        s.MarkerFaceAlpha = 0.4;
        hcb = colorbar('southoutside');
        hcb.Label.String= 'Likelihood L2';
        % title(sprintf('thingID = %d, zQSO = %.2f', selected_thing_ids, z_qso), 'FontSize', 20, 'Interpreter','latex');
        xlabel('$z_{CIV}$', 'FontSize', 20, 'Interpreter','latex');
        % xlabel('$\sigma$', 'FontSize', 20, 'Interpreter','latex');
        ylabel('$\log N_{CIV}$', 'FontSize', 20, 'Interpreter','latex');
        xlim([min(this_z_c4), max(this_z_c4)])
      
      
        
     
        
        fid = sprintf('EW/ID-%s.png',   ID{quasar_ind});
        % fid = sprintf('TN-Posterior-%s/TN-ID%s.pdf', training_set_name,  ID{quasar_ind});
        saveas(fig, fid, 'png');
        % fidpdf = sprintf('FP/%s.pdf',fid);
        % exportgraphics(fig, fid,'ContentType','vector')
    end
end
