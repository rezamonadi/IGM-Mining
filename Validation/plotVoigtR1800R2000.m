
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
    'model_posteriors', 'p_no_c4',...
    'p_L1', 'map_z_c4L1', 'map_N_c4L1', 'map_sigma_c4L1',...
    'map_z_c4L2', 'map_N_c4L2', 'map_sigma_c4L2'};

filename = sprintf('%s/processed_qsos_%s', ...
    processed_directory(release), ...
    training_set_name);
processed=matfile(filename);
load('data/C4_catalogs/Cooksey_C4_cat/processed/CIV-cat.mat','c4_QSO_ID','Z_c4','NCIV');

% load QSO model from training release
variables_to_load = {'rest_wavelengths', 'mu', 'M'};
load(sprintf('%s/learned_model',...
processed_directory(processed.training_release)), variables_to_load{:});

% load C4 samples from training release
variables_to_load = {'sigma_samples', 'offset_z_samples', 'log_nciv_samples', 'nciv_samples'};
load(sprintf('%s/civ_samples-%s', processed_directory(processed.training_release), training_set_name), ...
     variables_to_load{:});
     
% load preprocessed QSOs
variables_to_load = {'all_wavelengths', 'all_flux', 'all_noise_variance', ...
                     'all_pixel_mask'};
load(sprintf('%s/preloaded_qsos', processed_directory(release)), ...
     variables_to_load{:});

% building samples-> Z and sigma
% sample_sigma_c4 = min_sigma + (max_sigma-min_sigma)*offset_sigma_samples;




test_ind = processed.test_ind;
all_wavelengths    =    all_wavelengths(test_ind);
all_flux           =           all_flux(test_ind);
all_noise_variance = all_noise_variance(test_ind);
all_pixel_mask     =     all_pixel_mask(test_ind);
catalog = load(sprintf('%s/catalog', processed_directory(release)));
z_qsos = catalog.all_zqso(test_ind);
load(sprintf('%s/filter_flags', processed_directory(release)), ...
        'filter_flags');
p_L1 = processed.p_L1;
p_c4s = 1- processed.p_no_c4 -(1e-5)*p_L1;

sample_log_likelihoods_c4L2s = processed.sample_log_likelihoods_c4L2;
map_z_c4L2s = processed.map_z_c4L2;
map_N_c4L2s = processed.map_N_c4L2;
map_sigma_c4L2s = processed.map_sigma_c4L2;

map_z_c4L1s = processed.map_z_c4L1;
map_N_c4L1s = processed.map_N_c4L1;
map_sigma_c4L1s = processed.map_sigma_c4L1;

        % Testing Voigt profile 
dir = 'VoigtTest';
% dir = sprintf('TP-Posterior-N-sigma-%s',training_set_name);
mkdir(dir);
count =0;
for quasar_ind=1:10000
    
    %fprintf('quasar_ind:%d\n',quasar_ind);
    this_wavelengths    =    all_wavelengths{quasar_ind};
    this_wavelengths    =              this_wavelengths';
    this_flux           =           all_flux{quasar_ind}; 
    this_flux           =                     this_flux';
    this_noise_variance = all_noise_variance{quasar_ind};
    this_noise_variance =           this_noise_variance';
    this_pixel_mask     =     all_pixel_mask{quasar_ind};
    this_pixel_mask     =               this_pixel_mask';
    % convert to QSO rest frame
  
    this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qsos(quasar_ind));
        
    unmasked_ind = (this_rest_wavelengths >= min_lambda) & ...
        (this_rest_wavelengths <= max_lambda);
    % keep complete copy of equally spaced wavelengths for absorption
    % computation
    this_unmasked_wavelengths = this_wavelengths(unmasked_ind);

    % [mask_ind] remove flux pixels with pixel_mask; pixel_mask is defined
    % in read_spec_dr7.m
    ind = unmasked_ind & (~this_pixel_mask);

    this_wavelengths      =      this_wavelengths(ind);
    this_rest_wavelengths = this_rest_wavelengths(ind);
    this_flux             =             this_flux(ind);
    this_noise_variance   =   this_noise_variance(ind);

    % interpolate model onto given wavelengths
    mu_interpolator = ...
        griddedInterpolant(rest_wavelengths,        mu,        'linear');
    this_mu = mu_interpolator( this_rest_wavelengths);


    min_z_c4s = min_z_c4(this_wavelengths, z_qsos(quasar_ind));
    max_z_c4s = max_z_c4(this_wavelengths, z_qsos(quasar_ind));
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
    if(num_systems>0 & p_c4s(quasar_ind)<0.5) % FN
        % if(num_systems>0 & p_c4s(quasar_ind)>0.8 & filter_flags(quasar_ind)==0) % TP
        count=count+1;
        if count>20
            break
        end
        
        % fprintf('TP:%d\n',count);
        fprintf('FN:%d\n',count);
        fig = figure('visible', 'off', 'Position', [0,0,2024,1800]);
        clf();
        % subplot('position', [0.05 0.49 0.90 0.45]);
% construct dla_mu_map
        num_lines=2;
        absorptionL2R2000 = voigt(padded_wavelengths, map_z_c4L2s(quasar_ind), ...
            10^map_N_c4L2s(quasar_ind), num_lines, map_sigma_c4L2s(quasar_ind));
        
       

        absorptionL2R1800 = voigtR1800(padded_wavelengths, map_z_c4L2s(quasar_ind), ...
            10^map_N_c4L2s(quasar_ind), num_lines, map_sigma_c4L2s(quasar_ind));
        
       
        % Testing Voigt profile 
        
        all_c4_NCIV_test=all_c4_NCIV(test_ind);
        all_c4_Z_test=all_z_c4(test_ind);
        this_c4s = NCIV(this_systems); % cooksey's found C4s
        this_Zs  = Z_c4(this_systems); % Cooksey's found Zs
        
        matched_Z_c4L1= this_Zs(abs(this_Zs -map_z_c4L1s(quasar_ind))==min(abs(this_Zs-map_z_c4L1s(quasar_ind))));
        matched_N_c4L1= this_c4s(abs(this_Zs -map_z_c4L1s(quasar_ind))==min(abs(this_Zs-map_z_c4L1s(quasar_ind))));

        matched_Z_c4L2= this_Zs(abs(this_Zs -map_z_c4L2s(quasar_ind))==min(abs(this_Zs-map_z_c4L2s(quasar_ind))));
        matched_N_c4L2= this_c4s(abs(this_Zs -map_z_c4L2s(quasar_ind))==min(abs(this_Zs-map_z_c4L2s(quasar_ind))));
        
        absorptionL2R2000 = absorptionL2R2000(ind);
        c4_muL2R2000 = this_mu.* absorptionL2R2000;

        absorptionL2R1800 = absorptionL2R1800(ind);
        c4_muL2R1800 = this_mu.* absorptionL2R1800;
        this_z_c4 = (this_wavelengths / 1550) - 1;
        % xlim([]);
        % ylim([-0.1   5]);
        xlabel('(observed wavelengths $\lambda$ (\AA) / 1549 (\AA)) - 1', 'FontSize', 14, 'Interpreter','latex');
        ylabel('normalized flux $\mathbf{y}$',                            'FontSize', 14, 'Interpreter','latex');
        % plot(this_z_c4, c4_muL1R1800);
        % hold on
        plot(this_z_c4, c4_muL2R1800);
        hold on
        % plot(this_z_c4, c4_muL1R2000);
        % hold on
        plot(this_z_c4, c4_muL2R2000);
        % hold on
        % plot(this_z_c4, this_flux);
        legend( 'L2R1800', 'L2R2000')
        test_ind_c4 = all_z_c4>0;
        test_ind_c4= test_ind_c4(test_ind);
        
        % sub_title = sprintf('MAP(NL1):%.2f NL1C13:%.2f MAP(NL2):%.2f NL2C13:%.2f\nMAP(ZL1):%.2f ZL1C13:%.2f MAP(ZL2):%.2f ZL2C13:%.2f\nPL2:%.2f, PL1:%.2f\n',...
        %    map_N_c4L1s(quasar_ind), matched_N_c4L1, ...
        %    map_N_c4L2s(quasar_ind), matched_N_c4L2,...
        %    map_z_c4L1s(quasar_ind),matched_Z_c4L1,...
        %    map_z_c4L2s(quasar_ind),matched_Z_c4L2,...
        %    p_c4s(quasar_ind), p_L1(quasar_ind));
                  
        % [t,s] = title(sprintf('ID:%s Zqso:%.2f', ID{quasar_ind}, z_qsos(quasar_ind)), sub_title); 
        % t.FontSize=12;
        % s.FontSize=10;
        xlim([min(this_z_c4), max(this_z_c4)])
        xlabel('Z')
        ylabel('Normalized Flux')
        fid = sprintf('VoigtTest/%s.pdf', ID{quasar_ind});
        exportgraphics(fig, fid,'ContentType','vector')
    end
end