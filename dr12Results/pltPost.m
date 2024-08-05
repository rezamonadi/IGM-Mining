
fprintf('loading ... ')
filename = sprintf('%s/processed_qsos_dr12_N-1250-1610-S-35-115-nc-10k.mat', processed_directory(releaseTest));

load(filename)
N_CIV_DR12 =  savingCat.all_map_N_c4L2; 
Z_CIV_DR12 = savingCat.all_map_z_c4L2;
sigma_CIV_DR12 = savingCat.all_map_sigma_c4L2;
all_num_quasars = length(N_CIV_DR12)
p_all = savingCat.all_p_c4;

vs = {'N_CIV_DR12', 'Z_CIV_DR12', 'sigma_CIV_DR12', 'p_all' };
save('short_porseeed_dr12.mat',vs{:}, '-v7.3')


for tr= [0, 0.5, 0.65, 0.85, 0.95, 0.99]
    n = nnz(p_all>tr);
    fig = figure();
    h=histogram(reshape(N_CIV_DR12(p_all>tr), n,1), 'normalization', 'count')
    set(gca,'YScale','log')
    set(get(gca, 'XLabel'), 'String', '$N_{CIV} (cm^{-2})$', 'Interpreter','latex');
    ylabel('Count')
    set(gca, 'FontSize', 17)
    fprintf('sum(N): %d', sum(h.Values))
    exportgraphics(fig, sprintf('histPost/histDR12_NP%d.png',floor(tr*100)), 'Resolution', '800')

    fig = figure('visible', 'off');
    clf();
    h=histogram(reshape(Z_CIV_DR12(p_all>tr), n,1), 'normalization', 'count')
    set(gca,'YScale','log')
    set(get(gca, 'XLabel'), 'String', '$Z_{CIV}$', 'Interpreter', 'latex');
    ylabel('Count')
    set(gca, 'FontSize', 17)
    fprintf('sum(Z): %d', sum(h.Values))
    exportgraphics(fig, sprintf('histPost/histDR12_ZP%d.png',floor(tr*100)), 'Resolution', '800')

    fig = figure('visible', 'off');
    clf();
    h=histogram(reshape(sigma_CIV_DR12(p_all>tr), n, 1))
    set(gca,'YScale','log')
    set(get(gca, 'XLabel'), 'String', '$\sigma_{CIV} (cms^{-1})$', 'Interpreter','latex');
    set(gca, 'FontSize', 17)
    ylabel('Count')
    fprintf('sum(Sigma): %d', sum(h.Values))
    ylim([2000, max(h.Values)])
    exportgraphics(fig, sprintf('histPost/histDR12_SigmaP%d.png',floor(tr*100)), 'Resolution', '800')


end


h=histogram(reshape(Z_CIV_DR12(p_all>0.95), [],1), 'normalization', 'count')
h.BinWidth = 0.1;
x0 = h.BinEdges;
x0 = x0(2:end) - h.BinWidth/2;
y0 = h.Values;

% h=histogram(reshape(Z_CIV_DR12(p_all>0.85), [],1), 'normalization', 'count')
% h.BinWidth = 0.1;
% x1 = h.BinEdges;
% x1 = x1(2:end) - h.BinWidth/2;
% y1 = h.Values; 

% h=histogram(reshape(Z_CIV_DR12(p_all>0.65), [],1), 'normalization', 'count')
% h.BinWidth = 0.1;
% x2 = h.BinEdges;
% x2 = x2(2:end) - h.BinWidth/2;
% y2 = h.Values; 

z_qsos = all_zqso_dr12(test_ind);
h=histogram(z_qsos, 'normalization', 'count')
h.BinWidth = 0.1;
x3 = h.BinEdges;
x3 = x3(2:end) - h.BinWidth/2;
y3 = h.Values;


Z_C13 = reshape(all_z_civ_C13, [], 1);
Z_C13 = Z_C13(Z_C13>0);
h=histogram(Z_C13, 'normalization', 'count')
h.BinWidth = 0.1;
x4 = h.BinEdges;
x4 = x4(2:end) - h.BinWidth/2;
y4 = h.Values; 

fig = figure('visible', 'off');
clf();
plot(x0,y0, 'lineWidth', 2)
hold on


plot(x3,y3, 'lineWidth', 2, 'LineStyle', '--')
hold on

plot(x4,y4, 'lineWidth', 2, 'LineStyle', '-.')
hold on
set(gca,'YScale','log')
set(get(gca, 'XLabel'), 'String', '$Z$', 'Interpreter', 'latex');
ylabel('Count')
set(gca, 'FontSize', 17)
set(gca, 'YScale', 'log')
legend({'$Z_{CIV}(DR12)$',...
'$Z_{QSO}$', '$Z_{CIV}^{PM}$'}, 'interpreter','latex')
fprintf('sum(Z): %d', sum(h.Values))
exportgraphics(fig, 'histDR12_ZP95_curve_Z_c13_ZQSO.png', 'Resolution', '800')



% cleaning 
% clear N_CIV_DR12 Z_CIV_DR12 sigma_CIV_DR12


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % loading full 3D poster
p_c4       = savingCat.all_p_c4; 
p_c4L1     = savingCat.all_p_c4L1; 
p_no_c4     = savingCat.all_p_no_c4;
map_z_c4L2 = savingCat.all_map_z_c4L2;   
map_sigma_c4L2 = savingCat.all_map_sigma_c4L2;

%%% pTable%%%%%%%
numCIV=zeros(185425, 3);
j=0;
for tr=[0.65, 0.85, 0.95]
    j=j+1;
    for i=1:length(p_c4)
        numCIV(i, j)= nnz(p_c4(i,:)>=tr);
    end
end

for nAbs =1:7
    j=0;
    for tr=[0.65, 0.85, 0.95]
        j=j+1;
        
        fprintf('%d QSOs with %d Absorbers with P > %f\n', nnz(numCIV(:,j)==nAbs), nAbs, tr);
    end
    
end
fprintf('Sum(CIV 65)=%d\n',sum(numCIV(:,1)))
fprintf('Sum(CIV 85)=%d\n',sum(numCIV(:,2)))
fprintf('Sum(CIV 95)=%d\n',sum(numCIV(:,3)))
save('numCIV.mat', 'p_c4', 'numCIV', '-v7.3')

%%%   PC4 plot   %%%
fig = figure();
dat = p_c4(:,1);
dat = dat(~isnan(dat));
p = histogram(dat, 'Normalization', 'pdf');
x1 = p.BinEdges;
x1 = x1(2:end) - p.BinWidth/2;
y1 = p.Values; 

dat = p_c4(:,2);
dat = dat(~isnan(dat));
p = histogram(dat, 'Normalization', 'pdf');

x2 = p.BinEdges;
x2 = x2(2:end) - p.BinWidth/2;
y2 = p.Values; 
% hold on
dat = p_c4(:,3);
dat = dat(~isnan(dat));
p = histogram(dat, 'Normalization', 'pdf');

x3 = p.BinEdges;
x3 = x3(2:end) - p.BinWidth/2;
y3 = p.Values; 

dat = p_c4(:,4);
dat = dat(~isnan(dat));
p = histogram(dat, 'Normalization', 'pdf');
% p.FaceColor = 
% p.FaceAlpha = 0.6;
x4 = p.BinEdges;
x4 = x4(2:end) - p.BinWidth/2;
y4 = p.Values; 

clf();
l = plot(x1,y1); hold on
l.Color = [0.4660 0.6740 0.1880 0.7];
l.LineWidth = 2;

l=plot(x2,y2); hold on
l.LineWidth = 2;
l.Color = [0.6350 0.0780 0.1840 0.7];

l=plot(x3,y3); hold on
l.Color = [0 0.4470 0.7410 0.7];
l.LineWidth = 2;

l = plot(x4,y4); hold on
l.Color = [0.4940 0.1840 0.5560, 0.7];
l.LineWidth = 2;

xlabel('$P_n(M_D)$', 'Interpreter', 'latex')
legend({'1st search','2nd search', '3rd search', '4th search'},'Location', 'northwest')
ylabel('PDF')
set(gca, 'FontSize',17)
xlim([0,1])
exportgraphics(fig, 'pltPC4-v2/pDistAll.png', 'Resolution', 800);

% %-------------------------------------------------------
% fig = figure();
% p = histogram(p_c4(:,1));
% p.FaceColor = [0.4660 0.6740 0.1880];
% p.FaceAlpha = 0.7;
% xlabel('$P_n(Model_D)(1)$', 'Interpreter', 'latex')
% ylabel('Abundance')
% set(gca, 'FontSize',17)

% exportgraphics(fig, 'pltPC4-v2/pDist1.png', 'Resolution', 800);

% fig = figure();
% p = histogram(p_c4(:,2));
% p.FaceColor = [0.6350 0.0780 0.1840]; 
% p.FaceAlpha = 0.7;
% xlabel('$P_n(Model_D)(2)$', 'Interpreter', 'latex')
% ylabel('Abundance')
% set(gca, 'FontSize',17)

% exportgraphics(fig, 'pltPC4-v2/pDist2.png', 'Resolution', 800);

% fig = figure();
% p = histogram(p_c4(:,3));
% p.FaceColor = [0 0.4470 0.7410];
% p.FaceAlpha = 0.7;
% xlabel('$P_n(Model_D)(3)$', 'Interpreter', 'latex')
% ylabel('Abundance')
% set(gca, 'FontSize',17)

% exportgraphics(fig, 'pltPC4-v2/pDist3.png', 'Resolution', 800);

% fig = figure();
% p = histogram(p_c4(:,4));
% p.FaceColor = [0.4940 0.1840 0.5560];
% p.FaceAlpha = 0.7;
% xlabel('$P_n(Model_D)(4)$', 'Interpreter', 'latex')
% ylabel('Abundance')
% set(gca, 'FontSize',17)

% exportgraphics(fig, 'pltPC4-v2/pDist4.png', 'Resolution', 800);



%   PS plot   %%%
fig = figure();
dat = p_c4L1(:,1);
dat = dat(~isnan(dat));
p = histogram(dat);
x1 = p.BinEdges;
x1 = x1(2:end) - p.BinWidth/2;
y1 = p.Values; 

dat = p_c4L1(:,2);
dat = dat(~isnan(dat));
p = histogram(dat);

x2 = p.BinEdges;
x2 = x2(2:end) - p.BinWidth/2;
y2 = p.Values; 
% hold on
dat = p_c4L1(:,3);
dat = dat(~isnan(dat));
p = histogram(dat);

x3 = p.BinEdges;
x3 = x3(2:end) - p.BinWidth/2;
y3 = p.Values; 

dat = p_c4L1(:,4);
dat = dat(~isnan(dat));
p = histogram(dat);
% p.FaceColor = 
% p.FaceAlpha = 0.6;
x4 = p.BinEdges;
x4 = x4(2:end) - p.BinWidth/2;
y4 = p.Values; 

clf();
l = plot(x1,y1); hold on
l.Color = [0.4660 0.6740 0.1880 0.7];
l.LineWidth = 2;

l=plot(x2,y2); hold on
l.Color = [0.6350 0.0780 0.1840 0.7];
l.LineWidth = 2;

l=plot(x3,y3); hold on
l.LineWidth = 2;
l.Color = [0 0.4470 0.7410 0.7];

l = plot(x4,y4);
l.Color = [0.4940 0.1840 0.5560, 0.7];
l.LineWidth = 2;

xlabel('$P_n(Model_S)$', 'Interpreter', 'latex')
legend({'1st search','2nd search', '3rd search', '4th search'}, 'Location','northwest')
ylabel('Count')
set(gca, 'FontSize',17)
xlim([0,1])

exportgraphics(fig, 'pltPC4-v2/pL1DistAll.png', 'Resolution', 800);



%-------------------------------------------------------
% fig = figure();
% p = histogram(p_c4L1(:,1));
% p.FaceColor = [0.4660 0.6740 0.1880];
% p.FaceAlpha = 0.7;
% xlabel('$P_n(Model_S)(1)$', 'Interpreter', 'latex')
% ylabel('Abundance')
% set(gca, 'FontSize',17)

% exportgraphics(fig, 'pltPC4-v2/psDist1.png', 'Resolution', 800);

% fig = figure();
% p = histogram(p_c4L1(:,2));
% p.FaceColor = [0.6350 0.0780 0.1840]; 
% p.FaceAlpha = 0.7;
% xlabel('$P_n(Model_S)(2)$', 'Interpreter', 'latex')
% ylabel('Abundance')
% set(gca, 'FontSize',17)

% exportgraphics(fig, 'pltPC4-v2/psDist2.png', 'Resolution', 800);

% fig = figure();
% p = histogram(p_c4L1(:,3));
% p.FaceColor = [0 0.4470 0.7410];
% p.FaceAlpha = 0.7;
% xlabel('$P_n(Model_S)(3)$', 'Interpreter', 'latex')
% ylabel('Abundance')
% set(gca, 'FontSize',17)

% exportgraphics(fig, 'pltPC4-v2/psDist3.png', 'Resolution', 800);

% fig = figure();
% p = histogram(p_c4L1(:,4));
% p.FaceColor = [0.4940 0.1840 0.5560];
% p.FaceAlpha = 0.7;
% xlabel('$P_n(Model_S)(4)$', 'Interpreter', 'latex')
% ylabel('Abundance')
% set(gca, 'FontSize',17)

% exportgraphics(fig, 'pltPC4-v2/psDist4.png', 'Resolution', 800);



%   PNull plot   %%%
fig = figure();
dat = p_no_c4(:,1);
dat = dat(~isnan(dat));
p = histogram(dat);
x1 = p.BinEdges;
x1 = x1(2:end) - p.BinWidth/2;
y1 = p.Values; 

dat = p_no_c4(:,2);
dat = dat(~isnan(dat));
p = histogram(dat);

x2 = p.BinEdges;
x2 = x2(2:end) - p.BinWidth/2;
y2 = p.Values; 
% hold on
dat = p_no_c4(:,3);
dat = dat(~isnan(dat));
p = histogram(dat);

x3 = p.BinEdges;
x3 = x3(2:end) - p.BinWidth/2;
y3 = p.Values; 

dat = p_no_c4(:,4);
dat = dat(~isnan(dat));
p = histogram(dat);
% p.FaceColor = 
% p.FaceAlpha = 0.6;
x4 = p.BinEdges;
x4 = x4(2:end) - p.BinWidth/2;
y4 = p.Values; 

clf();
l = plot(x1,y1); hold on
l.Color = [0.4660 0.6740 0.1880 0.7];
l.LineWidth = 2;

l=plot(x2,y2); hold on
l.Color = [0.6350 0.0780 0.1840 0.7];
l.LineWidth = 2;

l=plot(x3,y3); hold on
l.Color = [0 0.4470 0.7410 0.7];
l.LineWidth = 2;

l = plot(x4,y4);
l.Color = [0.4940 0.1840 0.5560, 0.7];
l.LineWidth = 2;

xlabel('$P_n(Model_N)$', 'Interpreter', 'latex')
ylabel('$P_n(M_N)$', 'Interpreter', 'latex')
legend('1st search','2nd search', '3rd search', '4th search')
set(gca, 'FontSize',17)
exportgraphics(fig, 'pltPC4-v2/pNullDistAll.png', 'Resolution', 800);

%-------------------------------------------------------

%% Correltion between ambigous PC4s and other parameters   %%%

% Z(QSO)-PC4
z_qsos             =     all_zqso_dr12(test_ind);
fig = figure();
p_c4_1 = p_c4(:,1);
indPlt = (abs(p_c4_1 - 0.3) < 0.125);
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
mSN = nan(nnz(test_ind), max_civ);
allSN = nan(nnz(test_ind), max_civ);
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
    this_flux = this_flux(ind);
    this_z_1548 = (this_wavelengths / civ_1548_wavelength) - 1;
    this_z_1550 = (this_wavelengths / civ_1550_wavelength) - 1;
    for num_c4=1:7
        dz_Doppler = kms_to_z(sqrt(2)*5*map_sigma_c4L2(all_quasar_ind, num_c4)/1e5); % in km to z
        indIntegration = abs(this_z_1548 - map_z_c4L2(all_quasar_ind, num_c4))< dz_Doppler*(1+z_qso) | ...
                         abs(this_z_1550 - map_z_c4L2(all_quasar_ind, num_c4))< dz_Doppler*(1+z_qso);
        
        mSN(all_quasar_ind, num_c4) = mean((this_flux(indIntegration))./sqrt(thisNoise(indIntegration)));
        allSN(all_quasar_ind, num_c4) = mean(this_flux./sqrt(thisNoise));
    end
    
end

vs =  {'mSN', 'allSN'};
save('SN.mat', vs{:}, '-v7.3')
% SN(region)-PC4
fig = figure();
p_c4_1 = p_c4(:,1);
indPlt = (abs(p_c4_1 - 0.3) < 0.1);
h = histogram(mSN(:, 1), 'Normalization','pdf');
hold on 
% h.NumBins = 50;
% h.FaceAlpha = 0.6;
x1 = h.BinEdges;
x1 = x1(2:end) - h.BinWidth/2;
y1 = h.Values;
mSNpc41_30 = mSN(indPlt,1);
h = histogram(mSNpc41_30, 'Normalization', 'pdf');
x2 = h.BinEdges;
x2 = x2(2:end) - h.BinWidth/2;
y2 = h.Values;
hold on
h.NumBins = 50;
h.FaceAlpha = 0.6;
clf();
l = plot(x1,y1); hold on
l.Color = [0.4660 0.6740 0.1880 0.7];
l.LineWidth = 2;

l=plot(x2,y2); hold on
l.Color = [0.6350 0.0780 0.1840 0.7];
l.LineWidth = 2;

xlabel('$\langle S/N \rangle$', 'Interpreter', 'latex')
ylabel('PDF')
legend({'All QSOs', '$|P_1(M_D)-0.3|<0.1$'}, 'Interpreter', 'latex')
set(gca, 'FontSize',17)
xlim([0,20])
exportgraphics(fig, 'pltPC4-v2/hist_pc4_mSN_30.png', 'resolution', 800)


% SN(all)-PC4

% SN(region)-PC4
fig = figure();
p_c4_1 = p_c4(:,1);
indPlt = (abs(p_c4_1 - 0.3) < 0.1);
h = histogram(allSN(:, 1), 'Normalization','pdf');
hold on 
% h.NumBins = 50;
% h.FaceAlpha = 0.6;
x1 = h.BinEdges;
x1 = x1(2:end) - h.BinWidth/2;
y1 = h.Values;
allSNpc41_30 = allSN(indPlt,1);
h = histogram(allSNpc41_30, 'Normalization', 'pdf');
x2 = h.BinEdges;
x2 = x2(2:end) - h.BinWidth/2;
y2 = h.Values;
hold on
h.NumBins = 50;
h.FaceAlpha = 0.6;
clf();
l = plot(x1,y1); hold on
l.Color = [0.4660 0.6740 0.1880 0.7];
l.LineWidth = 2;
l=plot(x2,y2); hold on
l.Color = [0.6350 0.0780 0.1840 0.7];
l.LineWidth = 2;
xlabel('$S/N$', 'Interpreter', 'latex')
ylabel('PDF')
legend({'All QSOs', '$|P_1(M_D)-0.3|<0.1$'}, 'Interpreter', 'latex')
set(gca, 'YScale', 'log')
set(gca, 'FontSize',17)
xlim([0,20])
exportgraphics(fig, 'pltPC4-v2/hist_pc4_allSN_30.png', 'resolution', 800)



% scatter-SN-PC4
fig = figure();
x = reshape(p_c4(p_c4>0), [],1);
y = reshape(allSN(p_c4>0),[],1);
h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','off');

xlabel('S/N')
xlabel('$P_n(\mathcal{M}_D)$','Interpreter','latex')
set(gca, 'FontSize',17)
ylim([0 20])
exportgraphics(fig, 'pltPC4-v2/scatter_pc4_allSN.png', 'resolution', 800)

% mean allSN(pc4)
SNp = zeros(100,1);
ii=0;
for tr=0:0.01:0.99
    ii=ii+1;
    SNp(ii) = mean(mean(allSN((p_c4>tr) & (p_c4<tr+0.01))));
end

plot( 0.005:0.01:0.985, SNp(1:end-1))
xlabel('$P_n(\mathcal{M}_D)$','Interpreter','latex')
ylabel('SN')
set(gca, 'FontSize',17)
exportgraphics(fig, 'pltPC4-v2/SN_PC4.png', 'resolution', 800)

save('SN.mat', 'allSN', 'mSN', 'allSNpc41_30', 'mSNpc41_30', '-v7.3')

% ZCIV-PC4
fig = figure();
p_c4_1 = p_c4(:,1);
indPlt = (abs(p_c4_1 - 0.3) < 0.125);
x = p_c4_1(indPlt); y = map_z_c4L2(indPlt,1);
h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','off');
xlabel('$P_n(Model_{D}(1))$', 'Interpreter', 'latex')
ylabel('z_{CIV}')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/hist2D_pc41_zCIV.png', 'resolution', 800)

% ZCIV-zQSO-PC4
z_qsos             =     all_zqso_dr12(test_ind);
fig = figure();
p_c4_1 = p_c4(:,1);
indPlt = (abs(p_c4_1 - 0.3) < 0.125);
x = p_c4_1(indPlt); y = z_qsos(indPlt) - map_z_c4L2(indPlt,1);
h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','off');
xlabel('$P_n(Model_{D}(1))$', 'Interpreter', 'latex')
ylabel('z_{QSO} - z_{CIV}')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/hist2D_pc41_zCIV.png', 'resolution', 800)



% PS-PC4
fig = figure();
p_c4_1 = p_c4(:,1);
indPlt = (abs(p_c4_1 - 0.3) < 0.125);
x = p_c4_1(indPlt); y = p_c4L1(indPlt,1);
h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','off');
xlabel('$P_n(Model_{GP}(1))$', 'Interpreter', 'latex')
ylabel('$P_n(Model_S(1))$', 'Interpreter', 'latex')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/hist2D_pc41_pS.png', 'resolution', 800)


% Length-PC4
fig = figure();
p_c4_1 = p_c4(:,1);
indPlt = (abs(p_c4_1 - 0.48) < 0.125);
x = p_c4_1(indPlt); y = Length(indPlt)';
h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','on');
h.NumBins  = [20 50];
xlabel('$P_n(Model_{D}(1))$', 'Interpreter', 'latex')
ylabel('Length')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/hist2D_pc41_Length.png', 'resolution', 800)

% Large Sigma and SN
fig = figure();
x = reshape(map_sigma_c4L2, [], 1); y = reshape(allSN, [], 1);
p = reshape(p_c4, [], 1);
n = reshape(N_CIV_DR12, [], 1);
xp = x((p > 0.85) & (y<15) & (x>80e5)); 
yp = y((p > 0.85 )& (y<15) & (x>80e5));
h = histogram2(xp,yp,'DisplayStyle','tile','ShowEmptyBins','off');
% h.NumBins  = [50, 50];
xlabel('$\sigma$', 'Interpreter', 'latex')
ylabel('$\langle S/N \rangle$', 'Interpreter', 'latex')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/hist2D_SN-SigmaHighP85.png', 'resolution', 800)

% Large Sigma and SN
fig = figure();
x = reshape(map_sigma_c4L2, [], 1); y = reshape(allSN, [], 1);
p = reshape(p_c4, [], 1);
xp = x((p > 0.85) & (y < 15)); 
yp = y((p > 0.85) & (y < 15));
h = histogram2(xp,yp,'DisplayStyle','tile','ShowEmptyBins','off');
% h.NumBins  = [50, 50];
xlabel('$\sigma$', 'Interpreter', 'latex')
ylabel('$\langle S/N \rangle$', 'Interpreter', 'latex')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/hist2D_SN-SigmaAllP85.png', 'resolution', 800)


% mean S/N for all absorbers hist
fig = figure();
h = histogram(allSN(:,1));
h.NumBins =40;
xlim([0 20])
xlabel('mean(S/N)')
set(gca, 'FontSize',17)

exportgraphics(fig, 'pltPC4-v2/histmSN.png', 'resolution', 800)


fprintf('mean(SN(p>0.85):%.2f', nanmean(allSN(p_c4(:,1)>0.85, 1)));
fprintf('mean(SN(p>0.85 & N>15.5):%.2f', nanmean(allSN(p_c4(:,1)>0.85 & N_CIV_DR12(:,1)>15.5, 1)));
fprintf('mean(SN(p>0.85 & N>15.5 & s>100kms):%.2f', nanmean(allSN(p_c4(:,1)>0.85 & N_CIV_DR12(:,1)>15.5 & sigma_CIV_DR12(:,1)>100e5, 1)));

fprintf('med(SN(p>0.85):%.2f', nanmedian(allSN(p_c4(:,1)>0.85, 1)));
fprintf('med(SN(p>0.85 & N>15.5):%.2f', nanmedian(allSN(p_c4(:,1)>0.85 & N_CIV_DR12(:,1)>15.5, 1)));
fprintf('med(SN(p>0.85 & N>15.5 & s>100kms):%.2f', nanmedian(allSN(p_c4(:,1)>0.85 & N_CIV_DR12(:,1)>15.5 & sigma_CIV_DR12(:,1)>100e5, 1)));




