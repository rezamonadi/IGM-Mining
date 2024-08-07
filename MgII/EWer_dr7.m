filename ='data/dr7/processed/processed_qsos_tst_N-1250-1610-S-35-115-civWVL.mat';
fprintf('loading DR7 ...');
load(filename);
%%



num_quasars = nnz(test_ind);

all_wavelengths    =    all_wavelengths(test_ind);
all_flux           =           all_flux(test_ind);
all_pixel_mask     =     all_pixel_mask(test_ind);
all_sigma_pixel    =    all_sigma_pixel(test_ind);
z_qsos             =           all_zqso(test_ind);
num_quasars        =                numel(z_qsos);
all_noise_variance =   all_noise_variance(test_ind);
% % % preprocess model interpolants
% % % griddedInterpolant does an interpolation and gives a function handle
% % % based on the grided data. If the data is 2xD, {x,y} are the the same size as
% % % row  and columns. M is like M=f(x,y) and is like a matrix with each element
% % % M(i,j) = f(x(i), y(j))
mu_interpolator = ...
     griddedInterpolant(rest_wavelengths,        mu,        'linear');
M_interpolator = ...
     griddedInterpolant({rest_wavelengths, 1:k}, M,         'linear');

% % initialize results with nan
REW_1548_DR7_flux                        = nan(num_quasars, max_civ);
REW_1548_DR7_voigt                       = nan(num_quasars, max_civ);
ErrREW_1548_flux                         = nan(num_quasars, max_civ);
REW_1550_DR7_flux                        = nan(num_quasars, max_civ);
REW_1550_DR7_voigt                       = nan(num_quasars, max_civ);
ErrREW_1550_flux                         = nan(num_quasars, max_civ);
sigma_width = 5;    
parfor quasar_ind = 1:num_quasars
    tic;
    z_qso = z_qsos(quasar_ind);
    fprintf('processing quasar %i/%i (z_QSO = %0.4f) ...\n', ...
                              quasar_ind, num_quasars, z_qso);
    
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
        
        num_lines=1;
        absorptionL1_fine_1548= voigt_iP(padded_wavelengths_fine,...
            map_z_c4L2(quasar_ind,  num_c4), 10^map_N_c4L2(quasar_ind, num_c4),...
            num_lines, map_sigma_c4L2(quasar_ind, num_c4), padded_sigma_pixels_fine);
        absorptionL1_1548 = Averager(absorptionL1_fine_1548, nAVG, lenW_unmasked);
        absorptionL1_1548 = absorptionL1_1548(ind);

        absorptionL1_fine_1550 = voigt_iP_1(padded_wavelengths_fine,...
        map_z_c4L2(quasar_ind,  num_c4), 10^map_N_c4L2(quasar_ind, num_c4),...
        num_lines, map_sigma_c4L2(quasar_ind, num_c4), padded_sigma_pixels_fine);
        absorptionL1_1550 = Averager(absorptionL1_fine_1550, nAVG, lenW_unmasked);
        absorptionL1_1550 = absorptionL1_1550(ind);
        % Equivalent width calculation 
        % this_rest_wavelengths = this_rest_wavelengths(ind);
            
        % REW_1548
        dz_Doppler = kms_to_z(sqrt(2)*sigma_width*map_sigma_c4L2(quasar_ind, num_c4)/1e5); % in km to z
        dz_mid = (civ_1550_wavelength - civ_1548_wavelength)*0.5/civ_1548_wavelength; % cm/s
        indIntegration1 =  ...
                            (this_z_1548 < map_z_c4L2(quasar_ind, num_c4) + dz_Doppler*(1+z_qso)) & ...
                            (this_z_1548 > map_z_c4L2(quasar_ind, num_c4) - dz_Doppler*(1+z_qso)) & ...
                            (this_z_1548 < map_z_c4L2(quasar_ind, num_c4) + dz_mid*(1+z_qso));
    
        intThisFlux = this_flux(indIntegration1);
        intThisMu = this_mu(indIntegration1);

        if nnz(indIntegration1)>2
            REW_1548_DR7_flux(quasar_ind, num_c4) = trapz(emitted_wavelengths(this_wavelengths(indIntegration1),z_qso) ,...
                                                    1- intThisFlux./intThisMu);
            ErrREW_1548_flux(quasar_ind, num_c4)  = sqrt(sum(10^(2*pixel_spacing)*this_noise_variance(indIntegration1)./ ...
                                                        intThisMu.^2));  
        else
            fprintf('QSO-%d-c4-%d-npx:%d\n', quasar_ind, num_c4, nnz(indIntegration1));
        end

        REW_1548_DR7_voigt(quasar_ind, num_c4) = trapz(emitted_wavelengths(finer(this_unmasked_wavelengths, nAVG), z_qso),...
        1-absorptionL1_fine_1548);
        
        

        fprintf('QSO:%d, c4:%d, z_civ:%.2f, p(CIV):%.2f, REW_1548_flux:%.2f, REW_1548_voigt:%.2f\nerr:%.5f, sigma_width:%d\n',...
        quasar_ind, num_c4, map_z_c4L2(quasar_ind, num_c4), ...
        p_c4(quasar_ind, num_c4),  REW_1548_DR7_flux(quasar_ind, num_c4),...
        REW_1548_DR7_voigt(quasar_ind, num_c4), ErrREW_1548_flux(quasar_ind, num_c4), sigma_width);

        % REW_1550
        dz_Doppler = kms_to_z(sqrt(2)*sigma_width*map_sigma_c4L2(quasar_ind, num_c4)/1e5); % in km to z
        dz_mid = (civ_1550_wavelength - civ_1548_wavelength)*0.5/civ_1548_wavelength; % cm/s
        indIntegration2 = (this_z_1550 > map_z_c4L2(quasar_ind, num_c4) - dz_mid*(1+z_qso)) & ...
                            (this_z_1550 < map_z_c4L2(quasar_ind, num_c4) + dz_Doppler*(1+z_qso)) & ...
                            (this_z_1550 > map_z_c4L2(quasar_ind, num_c4) - dz_Doppler*(1+z_qso));
        intThisFlux = this_flux(indIntegration2);
        intThisMu = this_mu(indIntegration2);
        if (nnz(indIntegration2)>2)
            REW_1550_DR7_flux(quasar_ind, num_c4) = trapz(emitted_wavelengths(this_wavelengths(indIntegration2),z_qso) ,...
                                                    1- intThisFlux./intThisMu);
            ErrREW_1550_flux(quasar_ind, num_c4)  = sqrt(sum(10^(2*pixel_spacing)*this_noise_variance(indIntegration2)./ ...
                                                                                                    intThisMu.^2));
        else
            fprintf('QSO-%d-c4-%d-npx:%d\n', quasar_ind, num_c4, nnz(indIntegration2));

        end

        % fprintf('size(W):%d-%d,size(abs):%d-%d', size(this_wavelengths(indIntegration)),...
        % size(absorptionL1(indIntegration)));
        REW_1550_DR7_voigt(quasar_ind, num_c4) = trapz(emitted_wavelengths(finer(this_unmasked_wavelengths, nAVG), z_qso),...
        1-absorptionL1_fine_1550);
                                                                    
        fprintf('QSO:%d, c4:%d, z_civ:%.2f, p(CIV):%.2f, REW_1550_flux:%.2f, REW_1550_voigt:%.2f\nerr:%.5f, sigma_width:%d\n',...
        quasar_ind, num_c4, map_z_c4L2(quasar_ind, num_c4), ...
        p_c4(quasar_ind, num_c4),  REW_1550_DR7_flux(quasar_ind, num_c4),...
        REW_1550_DR7_voigt(quasar_ind, num_c4), ErrREW_1550_flux(quasar_ind, num_c4), sigma_width);




        % if quasar_ind<20
        %     fig = figure();
        %     p = stairs(this_z_1548, this_flux);
        %     p.Color = [0.1,0.1,0.8, 0.4];
        %     p.LineWidth=0.1;
        %     hold on
        %     p=plot(this_z_1548, this_mu);
        %     p.Color = [0.8,0.1,0.1, 0.4];
        %     p.LineWidth =0.1;
        %     hold on
        %     p=stairs(this_z_1548(indIntegration1), this_flux(indIntegration1));
        %     p.Color = [0, 0.1, 1];
        %     p.LineWidth = 2;
        %     hold on
            
        %     xline(map_z_c4L2(quasar_ind, num_c4))
        %     hold on 
        %     xline(map_z_c4L2(quasar_ind, num_c4) + dz_mid*(1+z_qso))
        %     exportgraphics(fig, sprintf('EW/pltEW1548/QSO-%d-c4-%d-sigma_width_%d.png',quasar_ind,...
        %                     num_c4, sigma_width), 'Resolution', 800);

        %     % fig = figure();
        %     % % p = stairs(this_rest_wavelengths, this_flux);
        %     % % p.Color = [0.1,0.1,0.8, 0.4];
        %     % % p.LineWidth=0.1;
        %     % % hold on
        %     % % p=plot(this_rest_wavelengths, this_mu);
        %     % % p.Color = [0.8,0.1,0.1, 0.4];
        %     % % p.LineWidth =0.1;
        %     % % hold on
        %     % p=stairs(this_rest_wavelengths(indIntegration2), this_flux(indIntegration2));
        %     % p.Color = [0, 0.1, 1];
        %     % hold on
        %     % p=plot(this_rest_wavelengths(indIntegration2), this_mu(indIntegration2));
        %     % p.Color = [1,0.1,0];
        %     % exportgraphics(fig, sprintf('EW/pltEW1550/QSO-%d-c4-%d-sigma_width_%d.png',quasar_ind,...
        %     %                 num_c4, sigma_width), 'Resolution', 800);
        % end
    % end

       

        
    end
    
end
vs = {'REW_1548_DR7_flux', 'REW_1548_DR7_voigt',...
'ErrREW_1548_flux', 'REW_1550_DR7_flux', 'REW_1550_DR7_voigt',...
'ErrREW_1550_flux'};
save(sprintf('EW/REW_DR7_sigma_width_%d.mat', sigma_width),vs{:} , '-v7.3');



