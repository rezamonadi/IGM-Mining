filename = sprintf('%s/processed_qsos_dr12_N-1250-1610-S-35-115-nc-10k.mat', processed_directory(releaseTest));
fprintf('loading DR12 ...');
% load(filename);

% p_c4                        = savingCat.all_p_c4;
% % sample_log_likelihoods_c4L2 = savingCat.all_sample_log_likelihoods_c4L2;
% map_z_c4L2                  = savingCat.all_map_z_c4L2;
% map_N_c4L2                  = savingCat.all_map_N_c4L2;
% map_sigma_c4L2              = savingCat.all_map_sigma_c4L2;   
% test_ind                    = savingCat.test_ind;
% vs = {'p_c4', 'map_z_c4L2', 'map_N_c4L2', 'map_sigma_c4L2', 'test_ind'}
% save('Wr/ShortProcessedDR12.mat', vs{:}, '-v7.3')
load('Wr/ShortProcessedDR12.mat')
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
REW_1550_DR12_voigt                        = nan(num_quasars, max_civ);
REW_1548_DR12_flux                        = nan(num_quasars, max_civ);
ErrREW_1548_flux                          = nan(num_quasars, max_civ);
    
for quasar_ind = 5
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
    for num_c4=1:7
        if isnan(p_c4(quasar_ind, num_c4))
            continue;
        end

        % Voigt REW 1548
        num_lines=1;
        absorptionL1_fine= voigt_iP(padded_wavelengths_fine,...
            map_z_c4L2(quasar_ind,  num_c4), 10^map_N_c4L2(quasar_ind, num_c4),...
            num_lines, map_sigma_c4L2(quasar_ind, num_c4), padded_sigma_pixels_fine);
        absorptionL1 = Averager(absorptionL1_fine, nAVG, lenW_unmasked);
        absorptionL1_1548 = absorptionL1(ind);
        absorptionL1_1548(isnan(absorptionL1_1548)) = 1;
        
        % Voigt REW 1550
        num_lines=1;
        absorptionL1_fine= voigt_iP_1(padded_wavelengths_fine,...
            map_z_c4L2(quasar_ind,  num_c4), 10^map_N_c4L2(quasar_ind, num_c4),...
            num_lines, map_sigma_c4L2(quasar_ind, num_c4), padded_sigma_pixels_fine);
        absorptionL1 = Averager(absorptionL1_fine, nAVG, lenW_unmasked);
        absorptionL1_1550 = absorptionL1(ind);
        absorptionL1_1550(isnan(absorptionL1_1550)) = 1;

        % Equivalent width calculation 
        % this_rest_wavelengths = this_rest_wavelengths(ind);
        dz_Doppler = kms_to_z(3*map_sigma_c4L2(quasar_ind, num_c4)/1e5); % in km to z
        dz_mid = (civ_1550_wavelength - civ_1548_wavelength)*0.5/civ_1548_wavelength; % cm/s
        % dz = min(dz_mid, dz_Doppler);
        % indIntegration = (this_z_1548< map_z_c4L2(quasar_ind, num_c4) + dz_mid) & ...
        indIntegration = (this_z_1548< map_z_c4L2(quasar_ind, num_c4) + dz_Doppler*(1+z_qso)) & ...
                         (this_z_1548> map_z_c4L2(quasar_ind, num_c4) - dz_Doppler*(1+z_qso));
        intThisFlux = this_flux(indIntegration);
        intThisMu   = this_mu(indIntegration);
        
        
        REW_1548_DR12_voigt(quasar_ind, num_c4) = trapz(emitted_wavelengths(this_wavelengths, z_qso),...
        1-absorptionL1_1548);
        REW_1550_DR12_voigt(quasar_ind, num_c4) = trapz(emitted_wavelengths(this_wavelengths, z_qso),...
        1-absorptionL1_1550);
        REW_1548_DR12_flux(quasar_ind, num_c4) = trapz(emitted_wavelengths(this_wavelengths(indIntegration),z_qso) ,...
        1- intThisFlux./intThisMu);
        ErrREW_1548_flux(quasar_ind, num_c4)  = 10^(pixel_spacing)*sqrt(sum(this_noise_variance(indIntegration)./ ...
                                                            intThisMu.^2));
        
        fprintf('QSO:%d, c4:%d, z_civ:%.2f, p(CIV):%.2f\nREW_1548_flux:%.2f(%.2f), REW_1548_voigt:%.2f\n\n',...
        quasar_ind, num_c4, map_z_c4L2(quasar_ind, num_c4), ...
        p_c4(quasar_ind, num_c4),  REW_1548_DR12_flux(quasar_ind, num_c4),ErrREW_1548_flux(quasar_ind, num_c4),...
        REW_1548_DR12_voigt(quasar_ind, num_c4));
        
    end
    
end
% vs = {'REW_1548_DR12_flux', 'REW_1548_DR12_voigt',...
% 'ErrREW_1548_flux1'};
% save('REW_1548_DR12.mat',vs{:} , '-v7.3');
% % DR =REW_1548_DR12_voigt./ REW_1550_DR12_voigt;




