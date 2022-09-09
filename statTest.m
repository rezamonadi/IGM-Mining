set_parameters_dr7
f = "processed_qsos_tst_6voigt-0-mask-1-prior-1-OccamRazor-1-nC4-30000-plt-0-MaskinP-0-fixedPriorCIV.mat";
% f = "processed_qsos_tst_6voigt-0-mask-1-prior-1-OccamRazor-1-nC4-30000-plt-0-MaskinP-1.mat";
% f = "processed_qsos_tst_6voigt-0-mask-1-prior-1-OccamRazor-1-nC4-30000-plt-0-MaskinP-5.mat";
filename = sprintf('%s/%s', processed_directory(release), f);
load(filename, 'p_c4', 'test_ind', 'map_z_c4L2', 'min_z_c4s', 'max_z_c4s', 'map_N_c4L2',...
        'map_sigma_c4L2');
num_quasars = nnz(test_ind);
% Saving Kathy's data                
zQSO_test = all_zqso(test_ind);

% comparing Z and N
Z_C13 = all_z_civ(test_ind,:);
Rate_test = all_RATING(test_ind, :);
dv= 350;
tr=0.85;



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
% EW                          = nan(num_quasars, max_civ);
    
% for quasar_ind = 1:num_quasars
% % for quasar_ind = 1:20
%     tic;
%     z_qso = z_qsos(quasar_ind);
%     fprintf('processing quasar %i/%i (z_QSO = %0.4f) ...', ...
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
%         [logspace(log10(min(this_unmasked_wavelengths)) - width * pixel_spacing, ...
%         log10(min(this_unmasked_wavelengths)) - pixel_spacing,...
%         width)';...
%         this_unmasked_wavelengths;...
%         logspace(log10(max(this_unmasked_wavelengths)) + pixel_spacing,...
%         log10(max(this_unmasked_wavelengths)) + width * pixel_spacing,...
%         width)'...
%         ];

%     ind = (~this_pixel_mask(unmasked_ind));
%     lenW_unmasked = length(this_unmasked_wavelengths);
%     ind_not_remove = true(size(this_flux));
%     for num_c4=1:max_civ
%         num_lines=2;
%         absorptionL2_fine= voigt_iP(finer(padded_wavelengths, nAVG),...
%             map_z_c4L2(quasar_ind,  num_c4), 10^map_N_c4L2(quasar_ind, num_c4),...
%             num_lines, map_sigma_c4L2(quasar_ind, num_c4), finer(this_sigma_pixel, nAVG));
%         absorptionL2 = Averager(absorptionL2_fine, nAVG, lenW_unmasked);
%         absorptionL2 = absorptionL2(ind);
%         % Equivalent width calculation 
%         EW(quasar_ind, num_c4) = trapz(this_wavelengths(ind), 1-absorptionL2);
%     end
%     fprintf(' took %0.3fs.\n', toc);
% end
% save('EW_DR7.mat', 'EW', '-v7.3');
load('EW_DR7.mat', 'EW');




for thisMaxCiv = 1:7
    PM=0;
    GP=0;
    PM1GP1=0; 
    PM0GP1 = 0;
    PM1GP0 = 0;
    PMn1GP1 = 0; 
    lw =5;
    dNdZ_PM0GP1 = nan([num_quasars, thisMaxCiv]);
    N_c4_PM1GP1 = nan([num_quasars, thisMaxCiv]);
    N_c4_PM0GP1 = nan([num_quasars, thisMaxCiv]);
    EW_PM1GP1_high_z = nan([num_quasars, thisMaxCiv]);
    EW_PM1GP1_low_z = nan([num_quasars, thisMaxCiv]);
    EW_PM0GP1_high_z = nan([num_quasars, thisMaxCiv]);
    EW_PM0GP1_low_z = nan([num_quasars, thisMaxCiv]);
    EW_PM1GP1 = nan([num_quasars, thisMaxCiv]);
    EW_PM0GP1 = nan([num_quasars, thisMaxCiv]);
    for quasar_ind=1:num_quasars
        for i=1:17
            if(Rate_test(quasar_ind,i)>=2)
                PM = PM+1;
            end
        end
        for j=1:thisMaxCiv
            if (p_c4(quasar_ind, j)>=tr)
                GP = GP+1;
            end
        end

        for i=1:17
            if(Rate_test(quasar_ind,i)>=2)
                for j=1:thisMaxCiv
                    DZ = abs(Z_C13(quasar_ind, i) - map_z_c4L2(quasar_ind, j, 1));
                    if p_c4(quasar_ind, j)>=tr & DZ<kms_to_z(dv)*(1+zQSO_test(quasar_ind))
                        PM1GP1 = PM1GP1 + 1;
                        N_c4_PM1GP1(quasar_ind, j) = map_N_c4L2(quasar_ind, j);
                        % if map_z_c4L2(quasar_ind, j)> 2.5
                            % EW_PM1GP1_high_z(quasar_ind, j) = EW(quasar_ind, j);
                        % else
    
                            EW_PM1GP1(quasar_ind, j) = EW(quasar_ind, j);
                        % end
                        % indPM1GP1_in_GP(quasar_ind,j) = 1;
                        % indPM1GP1(PM1GP1,2) = j;
                        break;
                    end
                end
            end
        end

       
       for j=1:thisMaxCiv
            if (p_c4(quasar_ind, j)>=tr) & (all(Rate_test(quasar_ind, :)<2))
                    PM0GP1 = PM0GP1 + 1;
                    dNdZ_PM0GP1(quasar_ind, j) =1/(max_z_c4s(quasar_ind) - min_z_c4s(quasar_ind));
                    N_c4_PM0GP1(quasar_ind, j) = map_N_c4L2(quasar_ind, j);
                    % if map_z_c4L2(quasar_ind, j)> 2.5
                        EW_PM0GP1(quasar_ind, j) = EW(quasar_ind, j);
                    % else

                        % EW_PM0GP1_low_z(quasar_ind, j) = EW(quasar_ind, j);
                    % end
            end
        end

    %    

    end
    
    % fig = figure('visible', 'off');
    % clf();
    % p = histogram(reshape(dNdZ_PM0GP1, num_quasars*thisMaxCiv, 1), 'Normalization', 'probability');
    % p.NumBins = 15;
    
    % exportgraphics(fig, sprintf('dNdZ-PM0GP1-max-%d.pdf', thisMaxCiv), 'Resolution', 800)

    % fig = figure('visible', 'off');
    % clf();
    % p = histogram(reshape(N_c4_PM0GP1, num_quasars*thisMaxCiv, 1), 'Normalization','probability');
    % p.NumBins = 25;
    % hold on
    % p = histogram(reshape(N_c4_PM1GP1, num_quasars*thisMaxCiv, 1), 'Normalization','probability');
    % p.NumBins = 25;

    % legend('N(PM0GP1)', 'N(PM1GP1)');
    % exportgraphics(fig, sprintf('N_c4_PM0GP1-max-%d.png', thisMaxCiv), 'Resolution', 800)

    % fig = figure('visible', 'off');
    % clf();
    % p = histogram(reshape(EW_PM0GP1_high_z, num_quasars*thisMaxCiv, 1), 'Normalization','probability');
    % p.NumBins = 25;
    % hold on
    % p = histogram(reshape(EW_PM1GP1_high_z, num_quasars*thisMaxCiv, 1), 'Normalization','probability');
    % p.NumBins = 25;

    % legend(sprintf('EW(PM0GP1), #%d', nnz(EW_PM0GP1_high_z)), sprintf('EW(PM1GP1), #%d', nnz(EW_PM0GP1_high_z)));

    % exportgraphics(fig, sprintf('EW_PM0GP1_high_z-max-%d.png', thisMaxCiv), 'Resolution', 800)


%     fig = figure('visible', 'off');
%     clf();
%     p1 = histogram(reshape(EW_PM0GP1, num_quasars*thisMaxCiv, 1));
%     p1.NumBins = 10;
%     hold on
   
%     p2 = histogram(reshape(EW_PM1GP1, num_quasars*thisMaxCiv, 1));
%     p2.NumBins = 10;
%     hold on 
%     p2.BinEdges = p1.BinEdges;
%     hold on 
%     y = p2.Values./(p1.Values + p2.Values);
%     h = p1.BinEdges;
%     H = (h(1:end-1) + h(2:end))*0.5;
 
%     legend(sprintf('EW(PM0GP1), #%d', PM0GP1), sprintf('EW(PM1GP1), #%d', PM1GP1));

%     exportgraphics(fig, sprintf('EW.png', thisMaxCiv), 'Resolution', 800)

%     fig = figure('visible', 'off');
%     clf();
%     stairs(H, y)
%     legend('PM1GP1/GP')
%     exportgraphics(fig, 'EW-Compeletness.png')
    fprintf('max:%d\nPM = %d\nGP = %d\nPM0GP1=%d\nPM1GP1 = %d\nPurity = %.2f, Completness = %.2f\n',...
            thisMaxCiv, PM, GP, PM0GP1, PM1GP1, 100*PM1GP1/GP, 100*PM1GP1/PM);
end

