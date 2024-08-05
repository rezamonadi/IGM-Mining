filename = sprintf('%s/processed_qsos_dr12_N-1250-1610-S-35-115-nc-10k.mat', processed_directory(releaseTest));
fprintf('loading DR12 ...');
load(filename);

p_c4                        = savingCat.all_p_c4;
% sample_log_likelihoods_c4L2 = savingCat.all_sample_log_likelihoods_c4L2;
map_z_c4L2                  = savingCat.all_map_z_c4L2;
map_N_c4L2                  = savingCat.all_map_N_c4L2;
map_sigma_c4L2              = savingCat.all_map_sigma_c4L2;   
test_ind                    = savingCat.test_ind;

num_quasars = nnz(test_ind);

all_wavelengths    =    all_wavelengths(test_ind);
all_flux           =           all_flux(test_ind);
all_pixel_mask     =     all_pixel_mask(test_ind);
all_sigma_pixel    =    all_sigma_pixel(test_ind);
z_qsos             =           all_zqso_dr12(test_ind);
num_quasars        =                numel(z_qsos);
all_noise_variance =   all_noise_variance(test_ind);
% % preprocess model interpolants
% % griddedInterpolant does an interpolation and gives a function handle
% % based on the grided data. If the data is 2xD, {x,y} are the the same size as
% % row  and columns. M is like M=f(x,y) and is like a matrix with each element
% % M(i,j) = f(x(i), y(j))
mu_interpolator = ...
    griddedInterpolant(rest_wavelengths,        mu,        'linear');
M_interpolator = ...
    griddedInterpolant({rest_wavelengths, 1:k}, M,         'linear');

% initialize results with nan
REW_1548_DR12_voigt   = nan(num_quasars, max_civ, 4);
REW_1548_DR12_flux    = nan(num_quasars, max_civ, 4);
ErrREW_1548_flux      = nan(num_quasars, max_civ, 4);
REW_1550_DR12_voigt   = nan(num_quasars, max_civ, 4);
REW_1550_DR12_flux    = nan(num_quasars, max_civ, 4);
ErrREW_1550_flux      = nan(num_quasars, max_civ, 4);

    
for quasar_ind = 1:num_quasars
    tic;
    z_qso = z_qsos(quasar_ind);
    fprintf('processing quasar %i/%i (z_QSO = %0.4f) ...\n', ...
                              quasar_ind, num_quasars, z_qso);
    
    this_wavelengths    =    all_wavelengths{quasar_ind};
    this_flux           =           all_flux{quasar_ind}; 
    this_pixel_mask     =     all_pixel_mask{quasar_ind};
    this_sigma_pixel    =     all_sigma_pixel{quasar_ind};  
    this_noise_variance = all_noise_variance{quasar_ind};

    % 
    % convert to QSO rest frame
    this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qso);
    
    unmasked_ind = (this_rest_wavelengths >= min_lambda) & ...
        (this_rest_wavelengths <= max_lambda);% & (this_sigma_pixel>0);
    % keep complete copy of equally spaced wavelengths for absorption
    % computation
    this_unmasked_wavelengths = this_wavelengths(unmasked_ind);
    this_unmasked_sigma_pixel = this_sigma_pixel(unmasked_ind); % avoiding mismathed sizes for padded variavles and NaNs
    % [mask_ind] remove flux pixels with pixel_mask; pixel_mask is defined
    % in read_spec_dr7.m
    ind = unmasked_ind & (~this_pixel_mask);
    this_wavelengths      =      this_wavelengths(ind);
    this_rest_wavelengths = this_rest_wavelengths(ind);
    this_flux             =             this_flux(ind);
    this_sigma_pixel      =      this_sigma_pixel(ind);
    this_noise_variance   = this_noise_variance(ind);
    % interpolate model onto given wavelengths
    this_mu = mu_interpolator( this_rest_wavelengths);
   
    
    
    padded_wavelengths_fine = ...
                            [logspace(log10(min(this_unmasked_wavelengths)) - width * pixel_spacing/(nAVG+1), ...
                            log10(min(this_unmasked_wavelengths)) - pixel_spacing/(nAVG+1),...
                            width)';...
                            finer(this_unmasked_wavelengths, nAVG)';...
                            logspace(log10(max(this_unmasked_wavelengths)) + pixel_spacing/(nAVG+1),...
                            log10(max(this_unmasked_wavelengths)) + width * pixel_spacing/(nAVG+1),...
                            width)'...
                            ];

    padded_sigma_pixels_fine = ...
                            [this_sigma_pixel(1)*ones(width,1);...
                            finer(this_sigma_pixel, nAVG)';...
                            this_sigma_pixel(end)*ones(width,1)];

    ind = (~this_pixel_mask(unmasked_ind));
    lenW_unmasked = length(this_unmasked_wavelengths);
    ind_not_remove = true(size(this_flux));
    this_z_1548 = (this_wavelengths / civ_1548_wavelength) - 1;
    this_z_1550 = (this_wavelengths / civ_1550_wavelength) - 1;
    for num_c4=1:max_civ
        if p_c4(quasar_ind, num_c4)<0.85
            break;
        end
       

        num_lines=1;
        absorptionL1_fine= voigt_iP(padded_wavelengths_fine,...
            map_z_c4L2(quasar_ind,  num_c4), 10^map_N_c4L2(quasar_ind, num_c4),...
            num_lines, map_sigma_c4L2(quasar_ind, num_c4), padded_sigma_pixels_fine);
        absorptionL1 = Averager(absorptionL1_fine, nAVG, lenW_unmasked);
        absorptionL1_1548 = absorptionL1(ind);

        absorptionL1_fine= voigt_iP_1(padded_wavelengths_fine,...
        map_z_c4L2(quasar_ind,  num_c4), 10^map_N_c4L2(quasar_ind, num_c4),...
        num_lines, map_sigma_c4L2(quasar_ind, num_c4), padded_sigma_pixels_fine);
        absorptionL1 = Averager(absorptionL1_fine, nAVG, lenW_unmasked);
        absorptionL1_1550 = absorptionL1(ind);
        % Equivalent width calculation 
        % this_rest_wavelengths = this_rest_wavelengths(ind);
        jj=0;
        for sigma_width = [1,2,3,4, 5]
            jj=jj+1;
             
            % REW_1548
            dz_Doppler = kms_to_z(sqrt(2)*sigma_width*map_sigma_c4L2(quasar_ind, num_c4)/1e5); % in km to z
            dz_mid = (civ_1550_wavelength - civ_1548_wavelength)*0.5/civ_1548_wavelength; % cm/s
            indIntegration = (this_z_1548 < map_z_c4L2(quasar_ind, num_c4) + dz_mid) & ...
                             (this_z_1548 < map_z_c4L2(quasar_ind, num_c4) + dz_Doppler*(1+z_qso)) & ...
                             (this_z_1548 > map_z_c4L2(quasar_ind, num_c4) - dz_Doppler*(1+z_qso));
            intThisFlux = this_flux(indIntegration);
            intThisMu = this_mu(indIntegration);
           
            REW_1548_DR12_voigt(quasar_ind, num_c4, jj) = trapz(emitted_wavelengths(this_wavelengths(indIntegration), z_qso),...
            1-absorptionL1_1548(indIntegration));
            REW_1548_DR12_flux(quasar_ind, num_c4, jj) = trapz(emitted_wavelengths(this_wavelengths(indIntegration),z_qso) ,...
            1- intThisFlux./intThisMu);
            ErrREW_1548_flux(quasar_ind, num_c4)  = sqrt(sum(10^(2*pixel_spacing)*this_noise_variance(indIntegration)./ ...
                                                            intThisMu.^2));                                                            
            fprintf('QSO:%d, c4:%d, z_civ:%.2f, p(CIV):%.2f, REW_1548_flux:%.2f, REW_1548_voigt:%.2f\nerr:%.5f, sigma_wdth:%d\n',...
            quasar_ind, num_c4, map_z_c4L2(quasar_ind, num_c4), ...
            p_c4(quasar_ind, num_c4),  REW_1548_DR12_flux(quasar_ind, num_c4, jj),...
            REW_1548_DR12_voigt(quasar_ind, num_c4, jj), ErrREW_1548_flux(quasar_ind, num_c4, jj), sigma_width);

            % REW_1550
            dz_Doppler = kms_to_z(sqrt(2)*sigma_width*map_sigma_c4L2(quasar_ind, num_c4)/1e5); % in km to z
            dz_mid = (civ_1550_wavelength - civ_1548_wavelength)*0.5/civ_1548_wavelength; % cm/s
            indIntegration = (this_z_1550 > map_z_c4L2(quasar_ind, num_c4) - dz_mid) & ...
                             (this_z_1550 < map_z_c4L2(quasar_ind, num_c4) + dz_Doppler*(1+z_qso)) & ...
                             (this_z_1550 > map_z_c4L2(quasar_ind, num_c4) - dz_Doppler*(1+z_qso));
            intThisFlux = this_flux(indIntegration);
            intThisMu = this_mu(indIntegration);
           
        
            % fprintf('size(W):%d-%d,size(abs):%d-%d', size(this_wavelengths(indIntegration)),...
            % size(absorptionL1(indIntegration)));
            REW_1550_DR12_voigt(quasar_ind, num_c4, jj) = trapz(emitted_wavelengths(this_wavelengths(indIntegration), z_qso),...
            1-absorptionL1(indIntegration));
            REW_1550_DR12_flux(quasar_ind, num_c4, jj) = trapz(emitted_wavelengths(this_wavelengths(indIntegration),z_qso) ,...
            1- intThisFlux./intThisMu);
            ErrREW_1550_flux(quasar_ind, num_c4)  = sqrt(sum(10^(2*pixel_spacing)*this_noise_variance(indIntegration)./ ...
                                                            intThisMu.^2));                                                            
            fprintf('QSO:%d, c4:%d, z_civ:%.2f, p(CIV):%.2f, REW_1550_flux:%.2f, REW_1550_voigt:%.2f\nerr:%.5f, sigma_wdth:%d\n',...
            quasar_ind, num_c4, map_z_c4L2(quasar_ind, num_c4), ...
            p_c4(quasar_ind, num_c4),  REW_1550_DR12_flux(quasar_ind, num_c4),...
            REW_1550_DR12_voigt(quasar_ind, num_c4), ErrREW_1550_flux(quasar_ind, num_c4), sigma_width);




            if quasar_ind<10
                fig = figure();
                p = stairs(this_rest_wavelengths, this_flux);
                p.Color = [0.1,0.1,0.8, 0.4];
                p.LineWidth=0.1;
                hold on
                p=plot(this_rest_wavelengths, this_mu);
                p.Color = [0.8,0.1,0.1, 0.4];
                p.LineWidth =0.1;
                hold on
                p=stairs(this_rest_wavelengths(indIntegration), this_flux(indIntegration));
                p.Color = [0, 0.1, 1];
                hold on
                p=plot(this_rest_wavelengths(indIntegration), this_mu(indIntegration));
                p.Color = [1,0.1,0];
                exportgraphics(fig, sprintf('EWdr12/pltEWQSO-v1/QSO-%d-c4-%d-SW_%d.png',quasar_ind,...
                              num_c4, sigma_width), 'Resolution', 800);
            end
        end
    end
    
end
vs = {'REW_1548_DR12_flux', 'REW_1548_DR12_voigt',...
'ErrREW_1548_flux', 'REW_1550_DR12_flux', 'REW_1550_DR12_voigt',...
'ErrREW_1550_flux'};
save('EW/REW_DR12.mat',vs{:} , '-v7.3');

% plotting 

load('EW/REW_1548_DR12.mat')

jj=0;
for sigma_width = [1, 2, 3, 5]
    jj=jj+1;

    rDew1 = 1- REW_1548_DR12_flux(:,:,jj)./REW_1548_DR12_voigt(:,:,jj);
    rDew_e = (REW_1548_DR12_voigt(:,:,jj) - REW_1548_DR12_flux(:,:,jj)) ./ ErrREW_1548_flux(:,:,jj);

    r1 = reshape(rDew1, [],1);
    r1 = r1( r1>quantile(r1, 0.001) & r1<quantile(r1,0.999));

    r3 = reshape(rDew_e, [],1);
    r3 = r3( r3>quantile(r3, 0.001) & r3<quantile(r3,0.999));

    fig = figure();
    histogram(r1, 'normalization', 'pdf');
    xlabel('$1 - \frac{REW(Flux)}{REW(Voigt)}$', 'interpreter','latex')
    ylabel('pdf')
    set(gca, 'FontSize',17)
    exportgraphics(fig, sprintf('EW/rREW_1548_voigt_flux-SW_%d.png', sigma_wdth), 'Resolution',800)

  
    fig = figure();
    histogram(r3, 'normalization', 'pdf');
    xlabel('$\frac{REW(Voigt) - REW(Flux)}{\sigma(REW(Flux))}$', 'interpreter','latex')
    ylabel('pdf')
    set(gca, 'FontSize',17)

    exportgraphics(fig, sprintf('EW/rREW_1548_voigt_flux_err_SW_%d.png', sigma_width), 'Resolution',800)



    fig = figure();
    rew1 = reshape(REW_1548_DR12_flux, [], 1);
    rew1 = rew1(rew1>quantile(rew1,.001) & rew1<quantile(rew1,0.999));
    p = histogram(rew1, 'normalization', 'pdf');
    p.FaceAlpha = 0.6;
    hold on 
    rew2 = reshape(REW_1548_DR12_voigt, [], 1);
    rew2 = rew2( rew2>quantile(rew2,0.001) & rew2<quantile(rew2,0.999));
    p = histogram(rew2, 'normalization', 'pdf');
    p.FaceAlpha = .6;
    hold on 
    xlabel('$REW_{1548A}$', 'interpreter', 'latex')
    ylabel('PDF')
    legend('Flux','Voigt')
    set(gca, 'FontSize',17)

    exportgraphics(fig, sprintf('EW/histREW_1548_flux_voigt_SW_%d.png',sigma_width), 'Resolution',800)


    fig = figure();
    rew1 = reshape(REW_1550_DR12_flux, [],1);
    rew1 = rew1(rew1>quantile(rew1,.001) & rew1<quantile(rew1,0.999));
    p = histogram(rew1, 'normalization', 'pdf');
    p.FaceAlpha = 0.6;
    hold on 
    rew2 = reshape(REW_1550_DR12_voigt, [],1);
    rew2 = rew2( rew2>quantile(rew2,0.001) & rew2<quantile(rew2,0.999));
    p = histogram(rew2, 'normalization', 'pdf');
    p.FaceAlpha = .6;
    hold on 
    xlabel('$REW_{1550A}$', 'interpreter', 'latex')
    ylabel('PDF')
    legend('Flux','Voigt')
    set(gca, 'FontSize',17)

    exportgraphics(fig, sprintf('EW/histREW_1550_flux_voigt_SW_%d.png', sigma_width), 'Resolution',800)


    fig = figure();
    clf();
    rew1 = reshape(REW_1548_DR12_voigt(p_all> 0.95), [],1);
    rew1 = rew1(rew1>quantile(rew1,.001) & rew1<quantile(rew1,0.999));


    rew2 = reshape(all_REW_1548_DR7(all_REW_1548_DR7 > 0), [], 1);
    rew2 = rew2(rew1>quantile(rew2,.001) & rew2<quantile(rew2,0.999));
    h = histogram(rew2, 'normalization', 'pdf');
    h.BinWidth = 0.1;
    x1 = x1(2:end) - h.BinWidth/2;
    y1 = h.Values;

    fig  = figure();
    clf();
    plot(x0, y0, 'LineWidth',2)
    hold on

    plot(x1, y1, 'LineWidth', 2)
    hold on

    xlabel('$REW_{1548A}$', 'interpreter', 'latex')
    ylabel('PDF')
    legend('PM(DR7)','GP(DR12)')
    set(gca, 'FontSize', 17)

    exportgraphics(fig, sprintf('EW/REW_1548_GP_PM_%d.png', sigma_width), 'Resolution',800)


    fig = figure();
    rew1 = reshape(REW_1550_DR12(p_all> 0.95), [],1);
    rew1 = rew1(rew1>quantile(rew1,.001) & rew1<quantile(rew1,0.999));


    rew2 = reshape(all_REW_1550_DR7(all_REW_1550_DR7 > 0), [], 1);
    rew2 = rew2(rew1>quantile(rew2,.001) & rew2<quantile(rew2,0.999));
    h = histogram(rew2, 'normalization', 'pdf');
    h.BinWidth = 0.1;
    x1 = x1(2:end) - h.BinWidth/2;
    y1 = h.Values;

    fig  = figure();
    clf();
    plot(x0, y0, 'LineWidth',2)
    hold on

    plot(x1, y1, 'LineWidth', 2)
    hold on

    xlabel('$REW_{1550A}$', 'interpreter', 'latex')
    ylabel('PDF')
    legend('PM(DR7)','GP(DR12)')
    set(gca, 'FontSize', 17)

    exportgraphics(fig, sprintf('EW/REW_1550_GP_PM_SW_%d.png', sigma_width) , 'Resolution',800)


    DR_1548_1550_DR12_voigt = REW_1548_DR12_voigt./REW_1550_DR12_voigt;
    h = histogram(reshape(DR_1548_1550_DR12_voigt(~isnan(DR_1548_1550_DR12_voigt)), [], 1), 'normalization', 'pdf');
    h.BinWidth = 0.1;
    xlabel('$REW_{1548}^{Voigt}/REW_{1550}^{Voigt}$', 'interpreter', 'latex')
    ylabel('PDF')

    exportgraphics(fig, sprintf('EW/DR_1550_1548_voigt.png', sigma_width), 'resolution', 800)



    DR_1548_1550_DR12_flux = REW_1548_DR12_flux./REW_1550_DR12_flux;

    h = histogram(reshape(DR_1548_1550_DR12_flux(~isnan(DR_1548_1550_DR12_flux)), [], 1), 'normalization', 'pdf');
    h.BinWidth = 0.1;
    xlabel('$REW_{1548}^{flux}/REW_{1550}^{flux}$', 'interpreter', 'latex')
    ylabel('PDF')

    exportgraphics(fig, sprintf('EW/DR_1550_1548_flux_SW_%d.png', sigma_width), 'resolution', 800)
end


