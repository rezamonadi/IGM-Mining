filename = sprintf('%s/processed_dr12.mat', processed_directory(releaseTest));
fprintf('loading DR12 ...');
load(filename);

p_c4                        = savingCat.all_p_c4;
sample_log_likelihoods_c4L2 = savingCat.all_sample_log_likelihoods_c4L2;
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
REW_1548_DR12_voigt                        = nan(num_quasars, max_civ);
REW_1548_DR12_flux                        = nan(num_quasars, max_civ);
ErrREW_1548_flux1                         = nan(num_quasars, max_civ);
ErrREW_1548_flux2                         = nan(num_quasars, max_civ);
    
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
   
    
    
    padded_wavelengths = ...
        [logspace(log10(min(this_unmasked_wavelengths)) - width * pixel_spacing, ...
        log10(min(this_unmasked_wavelengths)) - pixel_spacing,...
        width)';...
        this_unmasked_wavelengths;...
        logspace(log10(max(this_unmasked_wavelengths)) + pixel_spacing,...
        log10(max(this_unmasked_wavelengths)) + width * pixel_spacing,...
        width)'...
        ];

    ind = (~this_pixel_mask(unmasked_ind));
    lenW_unmasked = length(this_unmasked_wavelengths);
    ind_not_remove = true(size(this_flux));
    this_z_1548 = (this_wavelengths / civ_1548_wavelength) - 1;
    for num_c4=1:max_civ
        if p_c4(quasar_ind, num_c4)<0.85
            break;
        end
        num_lines=1;
        absorptionL1_fine= voigt_iP(finer(padded_wavelengths, nAVG),...
            map_z_c4L2(quasar_ind,  num_c4), 10^map_N_c4L2(quasar_ind, num_c4),...
            num_lines, map_sigma_c4L2(quasar_ind, num_c4), finer(this_sigma_pixel, nAVG));
        absorptionL1 = Averager(absorptionL1_fine, nAVG, lenW_unmasked);
        absorptionL1 = absorptionL1(ind);
        % Equivalent width calculation 
        % this_rest_wavelengths = this_rest_wavelengths(ind);
        dz_Doppler = kms_to_z(sqrt(2)*5*map_sigma_c4L2(quasar_ind, num_c4)/1e5); % in km to z
        dz_mid = (civ_1550_wavelength - civ_1548_wavelength)*0.5/civ_1548_wavelength; % cm/s
        % dz = min(dz_mid, dz_Doppler);
        % indIntegration = (this_z_1548< map_z_c4L2(quasar_ind, num_c4) + dz_mid) & ...
        indIntegration = (this_z_1548< map_z_c4L2(quasar_ind, num_c4) + dz_Doppler*(1+z_qso)) & ...
                         (this_z_1548> map_z_c4L2(quasar_ind, num_c4) - dz_Doppler*(1+z_qso));
        intThisFlux = this_flux(indIntegration);
        intThisMu   = this_mu(indIntegration);
        if nnz(indIntegration)>=3
            for ig=2:nnz(indIntegration)-2
                if (intThisFlux(ig-1)>intThisMu(ig-1)) & ...
                     (intThisFlux(ig)>intThisMu(ig)) & ...
                     (intThisFlux(ig+1)>intThisMu(ig+1)) 
                    intThisFlux(ig-1)=NaN;
                    intThisMu(ig-1)=NaN;
                end
            end
        else
            continue;
        end
        
        fprintf('size(W):%d-%d,size(abs):%d-%d', size(this_wavelengths(indIntegration)),...
        size(absorptionL1(indIntegration)));
        REW_1548_DR12_voigt(quasar_ind, num_c4) = trapz(emitted_wavelengths(this_wavelengths(indIntegration), z_qso),...
        1-absorptionL1(indIntegration));
        REW_1548_DR12_flux(quasar_ind, num_c4) = trapz(emitted_wavelengths(this_wavelengths(indIntegration),z_qso) ,...
        1- intThisFlux./intThisMu);
        ErrREW_1548_flux1(quasar_ind, num_c4)  = sqrt(sum(10^(2*pixel_spacing)*this_noise_variance(indIntegration)./ ...
                                                            intThisMu));
        ErrREW_1548_flux2(quasar_ind, num_c4)  = sqrt(sum(10^(2*pixel_spacing)*this_noise_variance(indIntegration)./ ...
                                                            intThisMu.^2));                                                            
        fprintf('QSO:%d, c4:%d, z_civ:%.2f, p(CIV):%.2f, REW_1548_flux:%.2f, REW_1548_voigt:%.2f\nerr1:%.5f, err2:%.5f\n',...
        quasar_ind, num_c4, map_z_c4L2(quasar_ind, num_c4), ...
        p_c4(quasar_ind, num_c4),  REW_1548_DR12_flux(quasar_ind, num_c4),...
        REW_1548_DR12_voigt(quasar_ind, num_c4), ErrREW_1548_flux1(quasar_ind, num_c4), ErrREW_1548_flux2(quasar_ind, num_c4));
        % if quasar_ind<100
        %     fig = figure();
        %     p = stairs(this_rest_wavelengths, this_flux);
        %     p.Color = [0.1,0.1,0.8, 0.4];
        %     p.LineWidth=0.1;
        %     hold on
        %     p=plot(this_rest_wavelengths, this_mu);
        %     p.Color = [0.8,0.1,0.1, 0.4];
        %     p.LineWidth =0.1;
        %     hold on
        %     p=stairs(this_rest_wavelengths(indIntegration), this_flux(indIntegration));
        %     p.Color = [0, 0.1, 1];
        %     hold on
        %     p=plot(this_rest_wavelengths(indIntegration), this_mu(indIntegration));
        %     p.Color = [1,0.1,0];
        %     exportgraphics(fig, sprintf('pltDR12EW/QSO-%d-c4-%d.png',quasar_ind, num_c4), 'Resolution', 800);

    end
    
end
vs = {'REW_1548_DR12_flux', 'REW_1548_DR12_voigt',...
'ErrREW_1548_flux1', 'ErrREW_1548_flux2'};
save('REW_1548_DR12.mat',vs{:} , '-v7.3');


