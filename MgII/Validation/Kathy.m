vCut = 3000;
max_civ = 4;
dv_mask = 700;
set_parameters_dr7
fprintf('Building catalogs ...\n')
variables_to_load= {'all_QSO_ID', 'all_zqso', 'c4_QSO_ID',...
'all_N_civ', 'all_z_civ', 'all_EW1', 'all_EW2'};
load(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_load{:});
filename = sprintf('%s/processed_qsos_tst_%s.mat', ...
processed_directory(release), ...
testing_set_name);
variables_to_load = {'test_ind',  'p_c4','map_z_c4L2',...
                    'map_sigma_c4L2', 'map_N_c4L2', 'EW',...
                    'log_likelihoods_no_c4', ''};
load(filename, variables_to_load{:});
load(sprintf('%s/preloaded_qsos_%s.mat', processed_directory(release), training_set_name));

variables_to_load = {'all_wavelengths'};
load(sprintf('%s/preloaded_qsos_%s.mat', processed_directory(release), training_set_name), ...
    variables_to_load{:});

test_flux = all_flux(test_ind);
test_noise2 = all_noise_variance(test_ind);
test_pixel_mask = all_pixel_mask(test_ind);

num_quasars = nnz(test_ind);
% Saving Kathy's data                
zQSO_test = all_zqso(test_ind);

% comparing Z and N
ID = all_QSO_ID(test_ind);
Z_C13 = all_z_civ(test_ind,:);
N_C13 = all_N_civ(test_ind,:);
% EW1_c13 = all_EW1(test_ind, :);
% EW2_c13 = all_EW2(test_ind, :);
all_wavelengths = all_wavelengths(test_ind);
% fig = figure(); 

dv= 300;
tr=0.95;
nPM=0;
nGP=0;
TP=0; 
FP = 0;
Missed =0;
for lw=[0.01, 0.1, 0.5, 1, 2, 3, 4, 5]
for quasar_ind=1:num_quasars
    this_wavelengths = all_wavelengths{quasar_ind};
    specSN = mean(test_flux{quasar_ind}./sqrt(test_noise2{quasar_ind}));
    for i=1:17
        if(Z_C13(quasar_ind,i)>0)
            nPM = nPM+1;
        end
    end
    for j=1:4
        if (p_c4(quasar_ind, j)>tr)
            nGP = nGP+1;
        end
    end
    for i=1:17
        this_TP =TP;
        if(Z_C13(quasar_ind,i)>0)
            for j=1:4
                DZ = abs(Z_C13(quasar_ind, i) - map_z_c4L2(quasar_ind,j));
                if (DZ<kms_to_z(dv)*(1+zQSO_test(quasar_ind))) &...
                    (p_c4(quasar_ind, j)>tr)
                    TP = TP+1;
                    N_TP(TP) = N_C13(quasar_ind, i);
                    Z_TP(TP) = Z_C13(quasar_ind, i);
                    Zqso_TP(TP) = zQSO_test(quasar_ind);
                    SN_TP(TP) = specSN;
                    sigma_TP(TP) = map_sigma_c4L2(quasar_ind, j);
                    EW_TP(TP) = EW(quasar_ind, j);
                    L_TP(TP) = log_likelihoods_no_c4(quasar_ind, j);
                    c4_pixel_ind = abs(this_wavelengths - (1+map_z_c4L2(quasar_ind, j))*1548.2040)<lw | ...
                    abs(this_wavelengths - (1+map_z_c4L2(quasar_ind, j))*1550.7810)<lw;
                    num_pixel_civ_TP(TP) = nnz(c4_pixel_ind);
                    
                end
            end
            if TP==this_TP
                Missed = Missed+1;
                N_Missed(Missed) = N_C13(quasar_ind, i);
                Z_Missed(Missed) = Z_C13(quasar_ind, i);
                Zqso_Missed(Missed) = zQSO_test(quasar_ind);
                SN_Missed(Missed) = specSN;
                sigma_Missed(Missed) = map_sigma_c4L2(quasar_ind, j);
                % EW_Missed(Missed) = EW1_c13(quasar_ind, i) + EW2_c13(quasar_ind, i);
                L_Missed(Missed) = mean(log_likelihoods_no_c4(quasar_ind, :));
                this_civ_pixel=0
                for j=1:4
                    c4_pixel_ind = abs(this_wavelengths - (1+map_z_c4L2(quasar_ind, j))*1548.2040)<lw | ...
                    abs(this_wavelengths - (1+map_z_c4L2(quasar_ind, j))*1550.7810)<lw;
                    this_civ_pixel = this_civ_pixel + nnz(c4_pixel_ind);
                end
                    num_pixel_civ_Missed(Missed) = this_civ_pixel/4;

            end
        end
    end
    for j=1:4
        if(p_c4(quasar_ind, j)>tr)
            DZ = abs(Z_C13(quasar_ind, :) - map_z_c4L2(quasar_ind, j));
            if all(DZ>kms_to_z(dv)*(1+zQSO_test(quasar_ind)))
                FP = FP+1;
                N_FP(FP) = map_N_c4L2(quasar_ind, j);
                Z_FP(FP) = map_z_c4L2(quasar_ind, j);
                Zqso_FP(FP) = zQSO_test(quasar_ind);
                SN_FP(FP) = specSN;
                sigma_FP(FP) = map_sigma_c4L2(quasar_ind, j);
                EW_FP(FP) = EW(quasar_ind, j);
                L_FP(FP) = log_likelihoods_no_c4(quasar_ind, j);
                c4_pixel_ind = abs(this_wavelengths - (1+map_z_c4L2(quasar_ind, j))*1548.2040)<lw | ...
                abs(this_wavelengths - (1+map_z_c4L2(quasar_ind, j))*1550.7810)<lw;
                num_pixel_civ_FP(FP) = nnz(c4_pixel_ind);

            end
        end
    end


end
Purity = TP/nGP;
Completeness= TP/nPM;
fprintf('Purity=%.2f\nCompletenss:%.2f\nnGP:%d\nnPM:%d\n', Purity, Completeness,... 
nPM, nGP);
end
% fig = figure();
% p=histogram(N_TP, 'Normalization', 'probability');
% p.NumBins=20;
% hold on
% p=histogram(N_FP, 'Normalization', 'probability');
% p.NumBins=20;
% p=histogram(N_Missed, 'Normalization', 'probability');
% p.NumBins=20;
% legend('TP','FP','Missed')
% xlabel('N')
% set(get(gca, 'YLabel'), 'String', 'PDF');

% exportgraphics(fig, sprintf('vCut%d/N.png', vCut), 'Resolution', 800)


% fig = figure();
% p=histogram(SN_TP, 'Normalization', 'probability');
% p.NumBins=20;
% hold on
% p=histogram(SN_FP, 'Normalization', 'probability');
% p.NumBins=20;
% p=histogram(SN_FP, 'Normalization', 'probability');
% p.NumBins=20;
% p=histogram(SN_Missed, 'Normalization', 'probability');
% p.NumBins=20;
% legend('TP','FP', 'Missed')
% xlabel('Signal/Noise')
% set(get(gca, 'YLabel'), 'String', 'PDF');

% exportgraphics(fig, sprintf('vCut%d/SN.png', vCut), 'Resolution', 800)


% fig = figure();
% p=histogram(sigma_TP, 'Normalization', 'probability');
% p.NumBins=20;
% hold on
% p=histogram(sigma_FP, 'Normalization', 'probability');
% p.NumBins=20;
% p=histogram(sigma_Missed, 'Normalization', 'probability');
% p.NumBins=20;
% legend('TP','FP','Missed')
% xlabel('Sigma')
% set(get(gca, 'YLabel'), 'String', 'PDF');

% exportgraphics(fig, sprintf('vCut%d/Sigma.png', vCut), 'Resolution', 800)


% fig = figure();
% p=histogram(Zqso_TP, 'Normalization', 'probability');
% p.NumBins=20;
% hold on
% p=histogram(Zqso_FP, 'Normalization', 'probability');
% p.NumBins=20;
% p=histogram(Zqso_Missed, 'Normalization', 'probability');
% p.NumBins=20;
% legend('TP','FP','Missed')
% xlabel('Zqso')
% set(get(gca, 'YLabel'), 'String', 'PDF');

% exportgraphics(fig, sprintf('vCut%d/zQSO.png', vCut), 'Resolution', 800)


% fig = figure();
% p=histogram(Z_TP, 'Normalization', 'probability');
% p.NumBins=20;
% hold on
% p=histogram(Z_FP, 'Normalization', 'probability');
% p.NumBins=20;
% p=histogram(Z_Missed, 'Normalization', 'probability');
% p.NumBins=20;
% legend('TP','FP', 'Missed')
% xlabel('Z_{CIV}')
% set(get(gca, 'YLabel'), 'String', 'PDF');

% exportgraphics(fig, sprintf('vCut%d/zCIV.png', vCut), 'Resolution', 800)


% fig = figure();
% p=histogram(EW_TP, 'Normalization', 'probability');
% p.NumBins=20;
% hold on
% p=histogram(EW_FP, 'Normalization', 'probability');
% p.NumBins=20;
% p=histogram(EW_missed, 'Normalization', 'probability');
% p.NumBins=20;
% legend('TP','FP', 'Missed')
% xlabel('EW_{CIV}')
% set(get(gca, 'YLabel'), 'String', 'PDF');
% exportgraphics(fig, sprintf('vCut%d/EW_CIV.png', vCut), 'Resolution', 800)


% fig = figure();
% p=histogram(L_TP, 'Normalization', 'probability');
% p.NumBins=20;
% hold on
% p=histogram(L_FP, 'Normalization', 'probability');
% p.NumBins=20;
% p=histogram(L_Missed, 'Normalization', 'probability');
% p.NumBins=20;
% legend('TP','FP', 'Missed')
% xlabel('log(L)')
% set(get(gca, 'YLabel'), 'String', 'PDF');

% exportgraphics(fig, sprintf('vCut%d/L.png', vCut), 'Resolution', 800)


fig = figure();
p=histogram(num_pixel_civ_TP, 'Normalization', 'probability');
p.NumBins=5;
hold on
p=histogram(num_pixel_civ_FP, 'Normalization', 'probability');
p.NumBins=5;
p=histogram(num_pixel_civ_Missed, 'Normalization', 'probability');
p.NumBins=5;
legend('TP','FP', 'Missed')
xlabel('numPixel')
set(get(gca, 'YLabel'), 'String', 'PDF');

exportgraphics(fig, sprintf('vCut%d/pixelCIV-lw%d.png', vCut, lw*100), 'Resolution', 800)



% fig = figure();
% x=0.001*speed_of_light*abs(Z_TP - Zqso_TP)./(1+Zqso_TP);
% min(x)
% p=histogram(x, 'Normalization', 'probability');
% p.NumBins=20;
% hold on
% x = 0.001*speed_of_light*abs(Z_FP - Zqso_FP)./(1+Zqso_FP);
% min(x)
% p=histogram(x, 'Normalization', 'probability');
% p.NumBins=20;
% x = 0.001*speed_of_light*abs(Z_Missed -Zqso_Missed)./(1+Zqso_Missed);
% min(x)
% p=histogram(x, 'Normalization', 'probability');
% p.NumBins=20;
% legend('TP','FP', 'Missed')
% % xlabel('Z_{CIV} - Z_{QSO}')
% xlabel('DV (km/s)')
% set(get(gca, 'YLabel'), 'String', 'PDF');
% exportgraphics(fig, sprintf('vCut%d/DV-qso-civ.png', vCut), 'Resolution', 800)
close all
