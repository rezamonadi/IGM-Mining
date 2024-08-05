filename = 'data/dr12/processed/processed_qsos_dr12_N-1250-1610-S-35-115-nc-10k.mat';
% filename = 'data/dr12/processed/processedDR12/processed_qsos_tst_DR12_S_1_E_2000.mat';
fprintf('loading DR12 ...\n');
load(filename);
% % loading full 3D poster
p_c4       = savingCat.all_p_c4; 
p_c4L1     = savingCat.all_p_c4L1; 
p_no_c4     = savingCat.all_p_no_c4;
map_z_c4L2 = savingCat.all_map_z_c4L2;   
map_sigma_c4L2 = savingCat.all_map_sigma_c4L2;
REW_1548_DR12 = savingCat.all_REW_1548_dr12;
p_all = savingCat.all_p_c4;

vs = {'p_all','REW_1548_DR12','all_REW_DR7'};
save('REW_all.mat', vs{:} , '-v7.3')

% plot REW(dr12) and REW(pm)

h=histogram(reshape(REW_1548_DR12(p_all> 0.95), [], 1), 'normalization', 'count')
h.BinWidth = 0.1;
x0 = h.BinEdges;
x0 = x0(2:end) - h.BinWidth/2;
y0 = h.Values;

h=histogram(reshape(all_REW_DR7(all_REW_DR7> 0),[],1), 'normalization', 'count')
h.BinWidth = 0.1;
x1 = h.BinEdges;
x1 = x1(2:end) - h.BinWidth/2;
y1 = h.Values;


fig = figure('visible', 'off');
clf();
plot(x0,y0, 'lineWidth', 2)
hold on

plot(x1,y1, 'lineWidth', 2)
hold on

set(gca,'YScale','log')
set(get(gca, 'XLabel'), 'String', '$REW_{1548A}$', 'Interpreter', 'latex');
ylabel('Count')
set(gca, 'FontSize', 17)
legend({'$REW_{1548A}(DR12, Voigt)$', '$REW_{1548A}(PM)$'}, 'interpreter','latex')
fprintf('sum(Z): %d', sum(h.Values))
exportgraphics(fig, 'histPost/histDR12voigt_DR7-REW.png', 'Resolution', '800')



%% pTable%%%%%%%
numCIV=nan(length(p_c4), 3);
j=0;
for tr=[0.65, 0.85, 0.95]
    j=j+1;
    for i=1:length(p_c4)
        numCIV(i, j)= nnz(p_c4(i,:)>=tr);
    end
end

numCIV_all=nan(8,3);
for nAbs =0:7
    j=0;
    for tr=[0.65, 0.85, 0.95]
        j=j+1;
        numCIV_all(nAbs+1, j) = nnz(numCIV(:,j)==nAbs);
        fprintf('%d QSOs with %d C4 with P > %f\n', nnz(numCIV(:,j)==nAbs), nAbs, tr);
    end
end

%% pTable%%%%%%%
numL1=nan(length(p_c4L1), 3);
j=0;
for tr=[0.65, 0.85, 0.95]
    j=j+1;
    for i=1:length(p_c4L1)
        numL1(i, j)= nnz(p_c4L1(i,:)>=tr);
    end
end

numL1_all=nan(8,3);
for nAbs =0:7
    j=0;
    for tr=[0.65, 0.85, 0.95]
        j=j+1;
        numL1_all(nAbs+1, j) = nnz(numL1(:,j)==nAbs);
        fprintf('%d QSOs with %d L1  with P > %f\n', nnz(numL1(:,j)==nAbs), nAbs, tr);
    end
end
fname = sprintf('%s/pTable.mat', processed_directory(releaseTest));
vs = {'numCIV_all', 'numL1_all'};
save(fname, vs{:}, '-v7.3')

%%%   PC4 plot   %%%

fig = figure();
p = histogram(p_c4(:,1));
p.FaceColor = [0.4660 0.6740 0.1880];
p.FaceAlpha = 0.6;
hold on

p = histogram(p_c4(:,2));
p.FaceColor = [0.6350 0.0780 0.1840]; 
p.FaceAlpha = 0.6;
hold on

p = histogram(p_c4(:,3));
p.FaceColor = [0 0.4470 0.7410];
p.FaceAlpha = 0.6;
hold on 

p = histogram(p_c4(:,4));
p.FaceColor = [0.4940 0.1840 0.5560];
p.FaceAlpha = 0.6;
xlabel('$P_n(Model_D)$', 'Interpreter', 'latex')
legend({'1st search','2nd search', '3rd search', '4th search'}, 'Location','northwest')
set(gca, 'FontSize',17)
xlim([0,1])
exportgraphics(fig, 'pltPC4-v2/pDistAll.png', 'Resolution', 800);

%-------------------------------------------------------
fig = figure();
p = histogram(p_c4(:,1));
p.FaceColor = [0.4660 0.6740 0.1880];
p.FaceAlpha = 0.7;
xlabel('$P_n(Model_D)(1)$', 'Interpreter', 'latex')
ylabel('Abundance')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/pDist1.png', 'Resolution', 800);

fig = figure();
p = histogram(p_c4(:,2));
p.FaceColor = [0.6350 0.0780 0.1840]; 
p.FaceAlpha = 0.7;
xlabel('$P_n(Model_D)(2)$', 'Interpreter', 'latex')
ylabel('Abundance')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/pDist2.png', 'Resolution', 800);

fig = figure();
p = histogram(p_c4(:,3));
p.FaceColor = [0 0.4470 0.7410];
p.FaceAlpha = 0.7;
xlabel('$P_n(Model_D)(3)$', 'Interpreter', 'latex')
ylabel('Abundance')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/pDist3.png', 'Resolution', 800);

fig = figure();
p = histogram(p_c4(:,4));
p.FaceColor = [0.4940 0.1840 0.5560];
p.FaceAlpha = 0.7;
xlabel('$P_n(Model_D)(4)$', 'Interpreter', 'latex')
ylabel('Abundance')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/pDist4.png', 'Resolution', 800);

%%%   PS plot   %%%

fig = figure();
p = histogram(p_c4L1(:,1));
p.FaceColor = [0.4660 0.6740 0.1880];
p.FaceAlpha = 0.6;
hold on

p = histogram(p_c4L1(:,2));
p.FaceColor = [0.6350 0.0780 0.1840]; 
p.FaceAlpha = 0.6;
hold on

p = histogram(p_c4L1(:,3));
p.FaceColor = [0 0.4470 0.7410];
p.FaceAlpha = 0.6;
hold on 

p = histogram(p_c4L1(:,4));
p.FaceColor = [0.4940 0.1840 0.5560];
p.FaceAlpha = 0.6;
xlabel('$P_n(Model_D)$', 'Interpreter', 'latex')
legend({'1st search','2nd search', '3rd search', '4th search'}, 'Location','northwest')
set(gca, 'FontSize',17)
xlim([0,1])
exportgraphics(fig, 'pltPC4-v2/psDistAll.png', 'Resolution', 800);

%-------------------------------------------------------
fig = figure();
p = histogram(p_c4L1(:,1));
p.FaceColor = [0.4660 0.6740 0.1880];
p.FaceAlpha = 0.7;
xlabel('$P_n(Model_S)(1)$', 'Interpreter', 'latex')
ylabel('Abundance')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/psDist1.png', 'Resolution', 800);

fig = figure();
p = histogram(p_c4L1(:,2));
p.FaceColor = [0.6350 0.0780 0.1840]; 
p.FaceAlpha = 0.7;
xlabel('$P_n(Model_S)(2)$', 'Interpreter', 'latex')
ylabel('Abundance')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/psDist2.png', 'Resolution', 800);

fig = figure();
p = histogram(p_c4L1(:,3));
p.FaceColor = [0 0.4470 0.7410];
p.FaceAlpha = 0.7;
xlabel('$P_n(Model_S)(3)$', 'Interpreter', 'latex')
ylabel('Abundance')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/psDist3.png', 'Resolution', 800);

fig = figure();
p = histogram(p_c4L1(:,4));
p.FaceColor = [0.4940 0.1840 0.5560];
p.FaceAlpha = 0.7;
xlabel('$P_n(Model_S)(4)$', 'Interpreter', 'latex')
ylabel('Abundance')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/psDist4.png', 'Resolution', 800);



%%%   PNull plot   %%%

fig = figure();
p = histogram(p_no_c4(:,1));
p.FaceColor = [0.4660 0.6740 0.1880];
p.FaceAlpha = 0.6;
hold on

p = histogram(p_no_c4(:,2));
p.FaceColor = [0.6350 0.0780 0.1840]; 
p.FaceAlpha = 0.6;
hold on

p = histogram(p_no_c4(:,3));
p.FaceColor = [0 0.4470 0.7410];
p.FaceAlpha = 0.6;
hold on 

p = histogram(p_no_c4(:,4));
p.FaceColor = [0.4940 0.1840 0.5560];
p.FaceAlpha = 0.6;
xlabel('$P_n(Model_{GP})$', 'Interpreter', 'latex')
legend('1st search','2nd search', '3rd search', '4th search')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/pgpDistAll.png', 'Resolution', 800);

%-------------------------------------------------------
min(p_no_c4)
max(p_no_c4)
fig = figure();
p = histogram(p_no_c4(:,1));
p.FaceColor = [0.4660 0.6740 0.1880];
p.FaceAlpha = 0.7;
xlabel('$P_n(Model_{GP}(1))$', 'Interpreter', 'latex')
ylabel('Abundance')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/pgpDist1.png', 'Resolution', 800);

fig = figure();
p = histogram(p_no_c4(:,2));
p.FaceColor = [0.6350 0.0780 0.1840]; 
p.FaceAlpha = 0.7;
p.NumBins=20;
xlabel('$P_n(Model_{GP}(2))$', 'Interpreter', 'latex')
ylabel('Abundance')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/pgpDist2.png', 'Resolution', 800);

fig = figure();
p = histogram(p_no_c4(:,3));
p.FaceColor = [0 0.4470 0.7410];
p.FaceAlpha = 0.7;
p.NumBins=20;

xlabel('$P_n(Model_{GP}(3))$', 'Interpreter', 'latex')
ylabel('Abundance')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/pgpDist3.png', 'Resolution', 800);

fig = figure();
p = histogram(p_no_c4(:,4));
p.FaceColor = [0.4940 0.1840 0.5560];
p.FaceAlpha = 0.7;
p.NumBins=20;
xlabel('$P_n(Model_{GP}(4))$', 'Interpreter', 'latex')
ylabel('Abundance')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/pgpDist4.png', 'Resolution', 800);


%%% Correltion between ambigous PC4s and other parameters   %%%

% Z(QSO)-PC4
z_qsos             =     all_zqso_dr12(test_ind);
fig = figure();
p_c4_1 = p_c4(:,1);
indPlt = (abs(p_c4_1 - 0.48) < 0.125);
x = p_c4_1(indPlt); y = z_qsos(indPlt);
h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','off');
xlabel('$P_n(Model_{D}(1))$', 'Interpreter', 'latex')
ylabel('z_{QSO}')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/hist2D_pc41_ZQSO.png', 'Resolution', 800)


% mean(S/N) for 350 km/s around MAP(z_CIV)
all_noise_variance = all_noise_variance(test_ind);
all_wavelengths = all_wavelengths(test_ind);
all_flux = all_flux(test_ind);
all_pixel_mask = all_pixel_mask(test_ind);
all_sigma_pixel = all_sigma_pixel(test_ind);
for all_quasar_ind=1:nnz(test_ind)
    z_qso = z_qsos(all_quasar_ind);
    this_wavelengths    =    all_wavelengths{all_quasar_ind};
    this_pixel_mask = all_pixel_mask{all_quasar_ind};
    this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qso);
    this_sigma_pixel   = all_sigma_pixel{all_quasar_ind};
    this_flux          = all_flux{all_quasar_ind};
    unmasked_ind = (this_rest_wavelengths >= min_lambda) & ...
        (this_rest_wavelengths <= max_lambda) & (this_sigma_pixel>0);
    % keep complete copy of equally spaced wavelengths for absorption
    % computation
    this_unmasked_wavelengths = this_wavelengths(unmasked_ind);
    
    % [mask_ind] remove flux pixels with pixel_mask; pixel_mask is defined
    % in read_spec_dr7.m
    ind = unmasked_ind & (~this_pixel_mask);
    this_wavelengths      =      this_wavelengths(ind);
    Length(all_quasar_ind) = nnz(ind);
    thisNoise = all_noise_variance{all_quasar_ind};
    thisNoise = thisNoise(ind);
    this_z_1548 = (this_wavelengths / 1548.2040) - 1;
    for num_c4=1:7
        dz_Doppler = kms_to_z(sqrt(2)*5*map_sigma_c4L2(all_quasar_ind, num_c4)/1e5); % in km to z
        indIntegration = abs(this_z_1548- map_z_c4L2(all_quasar_ind, num_c4))< dz_Doppler*(1+z_qso);
   
        SN(all_quasar_ind, num_c4) = mean(this_flux(indIntegration)./sqrt(thisNoise(indIntegration)));
    end
    
end
save('pltPC4-v2/mSN.mat', 'SN', '-v7.3')
% SN(region)-PC4
fig = figure();
p_c4_1 = p_c4(:,1);
indPlt = (abs(p_c4_1 - 0.33) < 0.2);
h = histogram(SN(:, 1));
hold on 
h.NumBins = 50;
h.FaceAlpha = 0.6;
SNpc50 = SN(indPlt,1);
h = histogram(SNpc50);
hold on
h.NumBins = 50;
h.FaceAlpha = 0.6;
xlabel('mean(S/N)')
ylabel('Abundance')
legend('All', '|P-0.33|<0.2')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/hist_pc41_mSN1.png', 'resolution', 800)
save('mSN.mat', 'SN', 'SNpc50', '-v7.3')


% % ZCIV-PC4
% z_qsos             =     all_zqso_dr12(test_ind);
% fig = figure();
% p_c4_1 = p_c4(:,1);
% indPlt = (abs(p_c4_1 - 0.48) < 0.125);
% x = p_c4_1(indPlt); y = map_z_c4L2(indPlt,1);
% h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','off');
% xlabel('$P_n(Model_{D}(1))$', 'Interpreter', 'latex')
% ylabel('z_{CIV}')
% set(gca, 'FontSize',17)

% exportgraphics(fig, 'pltPC4-v2/hist2D_pc41_zCIV.png', 'resolution', 800)



% % PS-PC4
% z_qsos             =     all_zqso_dr12(test_ind);
% fig = figure();
% p_c4_1 = p_c4(:,1);
% indPlt = (abs(p_c4_1 - 0.48) < 0.125);
% x = p_c4_1(indPlt); y = p_c4L1(indPlt,1);
% h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','off');
% xlabel('$P_n(Model_{GP}(1))$', 'Interpreter', 'latex')
% ylabel('$P_n(Model_S(1))$', 'Interpreter', 'latex')
% set(gca, 'FontSize',17)

% exportgraphics(fig, 'pltPC4-v2/hist2D_pc41_pS.png', 'resolution', 800)


% % % Length-PC4
% z_qsos             =     all_zqso_dr12(test_ind);
% fig = figure();
% p_c4_1 = p_c4(:,1);
% indPlt = (abs(p_c4_1 - 0.48) < 0.125);
% x = p_c4_1(indPlt); y = Length(indPlt)';
% h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','off');
% h.NumBins  = [20 50];
% xlabel('$P_n(Model_{D}(1))$', 'Interpreter', 'latex')
% ylabel('Length')
% set(gca, 'FontSize',17)

% exportgraphics(fig, 'pltPC4-v2/hist2D_pc41_Length.png', 'resolution', 800)


