
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
variables_to_load= {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7', ...
'all_QSO_ID', 'all_RA', 'all_DEC', 'all_zqso', 'EW1', 'EW2',...
'all_N_civ','all_z_civ', 'all_RATING', 'dla_QSO_ID','log_posteriors_dla',...
 'log_posteriors_no_dla' };
load(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_load{:});

% testing_set_name = sprintf('sigma-%d-%d-N-%d-%d', min_sigma/1e5, max_sigma/1e5, uniform_min_log_nciv*10 , uniform_max_log_nciv*10)
testing_set_name = 'sigma-15-65'

filename = sprintf('data/dr7/processed/processed_qsos_R%s.mat', testing_set_name);
% filename = '/home/reza/gpc/data/dr7/processed/processed_qsos_Rsigma-20-100-N-135-158-.mat'; 
% filename ='data/dr7/processed/processed_qsos_RRvoigt-10000Smp0-tr-80-20-vCutReza-5000.mat';
% filename ='data/dr7/processed/processed_qsos_Rprior-fixed.mat';
variables_to_load = {'training_release', 'training_set_name', ...
    'c4_catalog_name', 'prior_ind', 'release', ...
    'test_ind', 'prior_z_qso_increase', ...
    'max_z_cut', 'min_z_c4s', 'max_z_c4s', ...
    'log_priors_no_c4', 'log_priors_c4', ...
    'log_likelihoods_no_c4',  ...
    'sample_log_likelihoods_c4','log_likelihoods_c4'...
    'log_posteriors_no_c4', 'log_posteriors_c4', ...
    'model_posteriors', 'p_no_c4', 'p_c4', ...
    'map_z_c4', 'map_N_c4', 'map_sigma_c4', 'EqW1', 'EqW2', 'DoubletRatio'};

load(filename, variables_to_load{:});
load('data/C4_catalogs/Cooksey_C4_cat/processed/CIV-cat.mat','c4_QSO_ID','Z_c4','NCIV');

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
load(sprintf('%s/preloaded_qsos', processed_directory(release)), ...
     variables_to_load{:});

% building samples-> Z and sigma
% sample_sigma_c4 = min_sigma + (max_sigma-min_sigma)*offset_sigma_samples;
sample_sigma_c4 = sigma_samples;




test_ind = test_ind;
all_wavelengths    =    all_wavelengths(test_ind);
all_flux           =           all_flux(test_ind);
all_noise_variance = all_noise_variance(test_ind);
all_pixel_mask     =     all_pixel_mask(test_ind);
all_sigma_pixel    =    all_sigma_pixel(test_ind);

catalog = load(sprintf('%s/catalog', processed_directory(release)));
z_qsos = catalog.all_zqso(test_ind);
load(sprintf('%s/filter_flags', processed_directory(release)), ...
        'filter_flags');


p_c4s = p_c4;
sample_log_likelihoods_c4s = sample_log_likelihoods_c4;
map_N_c4s = map_N_c4;
map_z_c4s = map_z_c4;
map_sigma_c4s = map_sigma_c4;

% Testing Voigt profile 

dir = testing_set_name; 
mkdir(dir);
count =0;
nqso = nnz(test_ind)

for quasar_ind=1000:nqso
    % if ~mask_dz(quasar_ind)
    %     continue
    % end

    if count>100
       break
    end
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
    if(num_systems==0 & p_c4s(quasar_ind)>0.90) % FP
        % if(num_systems==0 & p_c4s(quasar_ind)<0.2 & filter_flags(quasar_ind)==0) % TN
        count=count+1;
        fprintf('FP:%d\n',count);
        fig = figure('visible', 'off');
        % fig = figure('visible', 'off', 'Position', [0,0,2024,1800]);
        clf();
        % subplot('position', [0.05 0.49 0.90 0.45]);
% construct dla_mu_map
        num_lines=2;
        % absorptionL2 = voigt(padded_wavelengths, map_z_c4L2s(quasar_ind), ...
        %     10^map_N_c4L2s(quasar_ind), num_lines, map_sigma_c4L2s(quasar_ind), ...
        %     this_sigma_pixel);
        % plot(1- absorptionL2)
       
        absorption = voigt(padded_wavelengths, map_z_c4s(quasar_ind), ...
        10^map_N_c4s(quasar_ind), 2, map_sigma_c4s(quasar_ind), this_sigma_pixel);
        
        % num_lines = 1;
        % absorptionL1 = voigt(padded_wavelengths, map_z_c4L1s(quasar_ind), ...
        %     10^map_N_c4L1s(quasar_ind), num_lines, map_sigma_c4L1s(quasar_ind), ...
        %     this_sigma_pixel);
        % Testing Voigt profile 
        
        % all_c4_NCIV_test=all_c4_NCIV(test_ind);
        % all_c4_Z_test=all_z_c4(test_ind);
        % this_c4s = NCIV(this_systems); % cooksey's found C4s
        % this_Zs  = Z_c4(this_systems); % Cooksey's found Zs
        % % s2 = scatter(this_Zs, this_c4s, 50, 'r', 'filled', 'DisplayName', 'C13');
        % % s2.MarkerFaceAlpha=1;
        % % s2.Marker= 'X';
        % % s2.MarkerEdgeColor='k';
        % % legend('GP', 'C13')
        % % hold on
        % matched_Z_c4L1= this_Zs(abs(this_Zs -map_z_c4L1s(quasar_ind))==min(abs(this_Zs-map_z_c4L1s(quasar_ind))));
        % matched_N_c4L1= this_c4s(abs(this_Zs -map_z_c4L1s(quasar_ind))==min(abs(this_Zs-map_z_c4L1s(quasar_ind))));

        % matched_Z_c4L2= this_Zs(abs(this_Zs -map_z_c4L2s(quasar_ind))==min(abs(this_Zs-map_z_c4L2s(quasar_ind))));
        % matched_N_c4L2= this_c4s(abs(this_Zs -map_z_c4L2s(quasar_ind))==min(abs(this_Zs-map_z_c4L2s(quasar_ind))));
        % % absorption = voigt(padded_wavelengths, this_Zs(1), ...
        %     % 10^this_c4s(1), num_lines, 1661700.9357883865);
            
        % absorptionL1 = absorptionL1(ind);
        absorption = absorption(ind);
        % c4_muL1 = this_mu.* absorptionL1;
        c4_mu = this_mu.* absorption;
               
        this_z_c4 = (this_wavelengths / 1550) - 1;

        % xlim([]);
        % subplot(2,1,1);
        xlabel('(observed wavelengths $\lambda$ (\AA) / 1549 (\AA)) - 1', 'FontSize', 14, 'Interpreter','latex');
        ylabel('normalized flux $\mathbf{y}$',                            'FontSize', 14, 'Interpreter','latex');
        
        plot(this_z_c4, c4_mu, 'Color', 'b');
        hold on
        
        TF = islocalmin(absorption, 'SamplePoints', this_unmasked_wavelengths);
        % Mins = absorptionL2(TF);
        % if(nnz(TF)>1)
        %     if (Mins(1)>Mins(2))
        %         1-Mins(1)>0.25*(1-Mins(2));
        %     end
        % else
        %     1-Mins(2)>0.25*(1-Mins(1));
        % end
        % fprintf('R:%.3f\n', (1-Mins(1))/(1-Mins(2)))
        % fprintf('Sep:%.4f\n', (max(this_unmasked_wavelengths(TF))-...
                % min(this_unmasked_wavelengths(TF))));
        % plot(this_z_c4(TF), absorptionL2(TF),'*', 'Color', 'r');
        % plot(this_z_c4, c4_muL2, 'Color', 'g');
        % hold on
        plot(this_z_c4, this_flux, 'Color', 'r');
        % hold on
        % legend('L2', 'mu*L2', 'Raw')
        % legend('L2',  'Raw')
        % legend( 'L1', 'L2', 'Raw')
        % plot( (1 + z_qsos(quasar_ind)) * (rest_wavelengths / 1549) - 1, mu, 'Color', 'g');
        % alpha(0.7)


        % test_ind_c4 = all_z_c4>0;
        % test_ind_c4= test_ind_c4(test_ind);
        % tit = sprintf('ID:%s\nMAP(N):%.3f,  N:%.2f, zqso=%.2f\nMAP(z):%.3f, Z:%.2f, sigma:%.2f, P:%.2f',...
        %             ID{quasar_ind}, map_N_c4L2s(quasar_ind),  matched_N_c4, z_qsos(quasar_ind), ...
        %             map_z_c4s(quasar_ind), matched_Z_c4, map_sigma_c4s(quasar_ind)/1e5, p_c4s(quasar_ind));
        % sub_title = sprintf('MAP(NL1):%.2f  MAP(NL2):%.2f\nMAP(ZL1):%.2f,    MAP(ZL2):%.2f\nPL2:%.2f, PL1:%.2f\n',...
        %    map_N_c4L1s(quasar_ind),...
        %    map_N_c4L2s(quasar_ind),...
        %    map_z_c4L1s(quasar_ind),...
        %    map_z_c4L2s(quasar_ind),...
        %    p_c4s(quasar_ind), p_L1(quasar_ind));
        % tit = sprintf('ID:%s\nN:%.2f, zqso=%.2f, Z:%.2f, sigma:%.2f, P:%.2f',...
        %             ID{quasar_ind}, this_c4s(1), z_qsos(quasar_ind), ...
        %             this_Zs(1), map_sigma_c4 p_c4(quasar_ind));
        tit = sprintf('ID:%s\nzqso=%.2f, EW(1548):%.2f, EW(1550):%.2f, DR:%.2f\nsigma:%.2fkm/s, N:%.2f, zCIV:%.2f, P:%.2f\n Dmin:%.4f',...
                     ID{quasar_ind}, z_qsos(quasar_ind), EqW1(quasar_ind,1),...
                     EqW2(quasar_ind,1), DoubletRatio(quasar_ind,1),...
                     map_sigma_c4(quasar_ind,1)/1e5, map_N_c4(quasar_ind,1),...
                     map_z_c4(quasar_ind,1), p_c4s(quasar_ind),...
                     (max(this_unmasked_wavelengths(TF))-...
                     min(this_unmasked_wavelengths(TF)))/(1+z_qsos(quasar_ind)));

        %  tit = sprintf('ID:%s\nzqso=%.2f, sigma:%.2fkm/s, N:%.2f\nzCIV:%.2f, P:%.2f\n Dmin:%.4f',...
        %              ID{quasar_ind}, z_qsos(quasar_ind), ...
        %              processed.map_sigma_c4L2(quasar_ind,1)/1e5, processed.map_N_c4L2(quasar_ind,1),...
        %              processed.map_z_c4L2(quasar_ind,1), p_c4s(quasar_ind),...
        %              (max(this_unmasked_wavelengths(TF))-...
        %              min(this_unmasked_wavelengths(TF))));                     
        title(tit)
        % [t,s] = title(sprintf('ID:%s Zqso:%.2f', ID{quasar_ind}, z_qsos(quasar_ind)), sub_title); 
        % t.FontSize=12;
        % s.FontSize=10;
        
        % text(1.7,2.5,tit)
        % ylim([-1 3]);
        % fid = num2str(quasar_ind);
        xlim([min(this_z_c4), max(this_z_c4)])
        % ylim([-0.1   1.1]);

        xlabel('Z')
        ylabel('Normalized Flux')
        % hold on 



        % subplot(2,1,2);
        % hold on
        % norm_sample_log_likelihoods = processed.sample_log_likelihoods_c4L2(quasar_ind, :);
        % norm_sample_log_likelihoods = norm_sample_log_likelihoods - max(norm_sample_log_likelihoods);
        % norm_sample_log_likelihoods = norm_sample_log_likelihoods - log(sum(exp(norm_sample_log_likelihoods)));

        % % s=scatter(sample_sigma_c4, log_nciv_samples,20,...
        %             %  norm_sample_log_likelihoods, 'filled', 'DisplayName', 'sample likelihoods');
        % s=scatter(sample_z_c4, log_nciv_samples,20,...
        %              norm_sample_log_likelihoods, 'filled', 'DisplayName', 'sample likelihoods');
        % s.MarkerFaceAlpha = 0.4;
        % hcb = colorbar('southoutside');
        % caxis([0.3*min(norm_sample_log_likelihoods) max(norm_sample_log_likelihoods)]);
        % hcb.Label.String= 'Likelihood L2';
        % % title(sprintf('thingID = %d, zQSO = %.2f', selected_thing_ids, z_qso), 'FontSize', 20, 'Interpreter','latex');
        % xlabel('$z_{CIV}$', 'FontSize', 20, 'Interpreter','latex');
        % % xlabel('$\sigma$', 'FontSize', 20, 'Interpreter','latex');
        % ylabel('$\log N_{CIV}$', 'FontSize', 20, 'Interpreter','latex');
        % xlim([min(this_z_c4), max(this_z_c4)])
        % % caxis([-100 0]);
        % hold on
        
        % % subplot(3,1,3);
        % % hold on
        % % norm_sample_log_likelihoods = processed.sample_log_likelihoods_c4L1(quasar_ind, :);
        % % norm_sample_log_likelihoods = norm_sample_log_likelihoods - max(norm_sample_log_likelihoods);
        % % norm_sample_log_likelihoods = norm_sample_log_likelihoods - log(sum(exp(norm_sample_log_likelihoods)));

        % % % s=scatter(sample_sigma_c4, log_nciv_samples,20,...
        % %             %  norm_sample_log_likelihoods, 'filled', 'DisplayName', 'sample likelihoods');
        % % s=scatter(sample_z_c4, log_nciv_samples,20,...
        % %             norm_sample_log_likelihoods, 'filled', 'DisplayName', 'sample likelihoods');
        % % s.MarkerFaceAlpha = 0.4;
        % % hcb = colorbar('southoutside');
        % % hcb.Label.String= 'Likelihood L1';
        % % % title(sprintf('thingID = %d, zQSO = %.2f', selected_thing_ids, z_qso), 'FontSize', 20, 'Interpreter','latex');
        % % xlabel('$z_{CIV}$', 'FontSize', 20, 'Interpreter','latex');
        % % % xlabel('$\sigma$', 'FontSize', 20, 'Interpreter','latex');
        % % ylabel('$\log N_{CIV}$', 'FontSize', 20, 'Interpreter','latex');
        % % xlim([min(this_z_c4), max(this_z_c4)])

        % % % caxis([-100 0]);
     
        
        fid = sprintf('%s/FP-ID-%s.png',  dir,  ID{quasar_ind});
        % fid = sprintf('TN-Posterior-%s/TN-ID%s.pdf', training_set_name,  ID{quasar_ind});
        % saveas(fig, fid, 'png');
        % fidpdf = sprintf('FP/%s.pdf',fid);
        % exportgraphics(fig, fid,'ContentType','vector')
        exportgraphics(fig, fid,'Resolution',400)
    end
end
