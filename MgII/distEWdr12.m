load('Wr/ShortProcessedDR12.mat')
sigma_width=4;
minpx = 1;


for sn=0
    % plotting 
    
    load(sprintf('Wr/Wr_DR12_sigma_width_%d_minpx_%d.mat', sigma_width, minpx))

%     % scatter Diff_REW_flux_voigt vs. W_r,1548
%     dEW_voigt_flux_dr12_error = (REW_1548_DR12_voigt - REW_1548_DR12_flux) ./ ErrREW_1548_flux;

%     y = reshape(dEW_voigt_flux_dr12_error, [],1);
%     x = reshape(REW_1548_DR12_flux, [], 1);

%     fig = figure();
%     s = scatter(x,y, 10, 'filled');
%     s.MarkerFaceAlpha = 0.4;
%     xlabel('${\rm W}_{r,1548}^{\rm GP, flux}$', 'interpreter','latex')
%     % ylabel('\frac{{\rm W}_{r,1548}^{\rm GP, flux} - {\rm W}_{r,1548}^{\rm GP, Voigt}}{\rm err}', 'interpreter', 'latex')
%     ylabel('(WrFlux-WrVoigt)/err')
%     exportgraphics(fig, sprintf('Wr/scatter_dREW_1548_voigt_flux-SW_%d_minpx_%d.png', sigma_width, minpx), 'Resolution',800)
%     stop

%     % scatter DR vs. P_C4
    DR_1548_1550_DR12_flux = REW_1548_DR12_flux./REW_1550_DR12_flux;
%     fig = figure();
%     y = reshape(DR_1548_1550_DR12_flux(~isnan(DR_1548_1550_DR12_flux)), [], 1);
%     x = reshape(p_c4(~isnan(DR_1548_1550_DR12_flux)), [], 1);
%     s =scatter(x,y, 10, 'filled');
%     s.MarkerFaceAlpha = 0.4;
%     ylabel('${\rm W}_{r,1548}^{\rm GP,flux}/{\rm W}_{r,1550}^{\rm GP, flux}$', 'interpreter', 'latex')
%     xlabel('P(M$_D$)', 'interpreter','latex')
%     exportgraphics(fig, sprintf('Wr/scatter_DR_flux_pc4_SW_%d_minpx_%d.png', sigma_width, minpx), 'resolution', 800)

%     % scatter DR vs. minpx
%     fig = figure();
%     y = reshape(DR_1548_1550_DR12_flux(~isnan(DR_1548_1550_DR12_flux)), [], 1);
%     x = reshape(REW_1548_DR12_flux_indIntegration(~isnan(DR_1548_1550_DR12_flux)), [], 1);
%     s =scatter(x,y, 10, 'filled');
%     s.MarkerFaceAlpha = 0.4;
%     ylabel('${\rm W}_{r,1548}^{\rm GP,flux}/{\rm W}_{r,1550}^{\rm GP, flux}$', 'interpreter', 'latex')
%     xlabel('min(pixel1)')
%     exportgraphics(fig, sprintf('Wr/scatter_DR_flux_minpx_SW_%d_minpx_%d.png', sigma_width, minpx), 'resolution', 800)

%     % scatter DR vs. W_r,1548
%     fig = figure();
%     y = reshape(DR_1548_1550_DR12_flux(~isnan(DR_1548_1550_DR12_flux)), [], 1);
%     x = reshape(REW_1548_DR12_flux(~isnan(DR_1548_1550_DR12_flux)), [], 1);
%     s =scatter(x,y, 10, 'filled');
%     s.MarkerFaceAlpha = 0.4;
%     ylabel('${\rm W}_{r,1548}^{\rm GP,flux}/{\rm W}_{r,1550}^{\rm GP, flux}$', 'interpreter', 'latex')
%     xlabel('${\rm W}_{r,1548}^{\rm GP,flux}$','interpreter','latex')
%     exportgraphics(fig, sprintf('Wr/scatter_DR_flux_Wr1548_SW_%d_minpx_%d.png', sigma_width, minpx), 'resolution', 800)


%     % scatter DR vs. <SNR>
%     load('SN.mat');
%     fig = figure();
%     y = reshape(DR_1548_1550_DR12_flux(~isnan(DR_1548_1550_DR12_flux)), [], 1);
%     x = reshape(mSN(~isnan(DR_1548_1550_DR12_flux)), [], 1);
%     s =scatter(x,y, 10, 'filled');
%     s.MarkerFaceAlpha = 0.4;
%     ylabel('${\rm W}_{r,1548}^{\rm GP,flux}/{\rm W}_{r,1550}^{\rm GP, flux}$', 'interpreter', 'latex')
%     xlabel('$\langle{\rm SNR}\rangle$','interpreter','latex')
%     exportgraphics(fig, sprintf('Wr/scatter_DR_mSNR_SW_%d_minpx_%d.png', sigma_width, minpx), 'resolution', 800)

%    % scatter DR vs. err
%     fig = figure();
%     y = reshape(DR_1548_1550_DR12_flux(~isnan(DR_1548_1550_DR12_flux)), [], 1);
%     x = reshape(ErrREW_1548_flux(~isnan(DR_1548_1550_DR12_flux)), [], 1);
%     s =scatter(x,y, 10, 'filled');
%     s.MarkerFaceAlpha = 0.4;
%     ylabel('${\rm W}_{r,1548}^{\rm GP,flux}/{\rm W}_{r,1550}^{\rm GP, flux}$', 'interpreter', 'latex')
%     xlabel('err(W$_{r,1548}$)', 'interpreter','latex')
%     exportgraphics(fig, sprintf('Wr/scatter_DR_errWr1548_SW_%d_minpx_%d.png', sigma_width, minpx), 'resolution', 800)

%     % % P-S/N hist2d of 
    load('SN.mat')
%     % fig = figure();
%     % x = reshape(REW_1550_DR12_voigt(p_c4 > 0), [],1);
%     % y = mSN;
%     % y = y(p_c4>0);
%     % y = reshape(y, [], 1);
%     % h = histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','off');
%     % h.NumBins =[100,100];
%     % xlabel('REW(1548A)')
%     % ylabel('$\langle S/N \rangle$', 'interpreter', 'latex')
%     % set(gca, 'FontSize',15)

%     % exportgraphics(fig, 'EW/hist2D_REW_mSNP85.png', 'Resolution', 800)

%     % % Sigma-S/N hist2d of 
%     % load('SN.mat')
%     % x1 = reshape(mSN(map_sigma_c4L2< 80e5 & mSN>0), [],1);
%     % x2 = reshape(mSN(map_sigma_c4L2>= 80e5 & map_sigma_c4L2 < 100e5 & mSN>0), [],1);
%     % x3 = reshape(mSN(map_sigma_c4L2>= 100e5 &  mSN>0), [],1);

%     % h = histogram(x1, 'normalization', 'pdf')
%     % h.BinWidth = 0.5;
%     % x1 = h.BinEdges;
%     % x1 = x1(2:end) - h.BinWidth/2;
%     % y1 = h.Values;


%     % h = histogram(x2, 'normalization', 'pdf')
%     % h.BinWidth = 0.5;
%     % x2 = h.BinEdges;
%     % x2 = x2(2:end) - h.BinWidth/2;
%     % y2 = h.Values;



%     % h = histogram(x3, 'normalization', 'pdf')
%     % h.BinWidth = 0.5;
%     % x3 = h.BinEdges;
%     % x3 = x3(2:end) - h.BinWidth/2;
%     % y3 = h.Values;

%     % fig = figure();

%     % p = plot(x1, y1)
%     % p.LineWidth = 2;
%     % hold on 
%     % plot(x2,y2, 'LineWidth', 2)
%     % hold on 
%     % plot(x3,y3, 'LineWidth', 2)

%     % legend({'$\sigma_{CIV}>80kms^{-1}$',...
%     %          '$80kms^{-1}\le\sigma_{CIV}<100kms^{-1}$',...
%     %          '$\sigma_{CIV}>100kms^{-1}$'}, 'interpreter', 'latex')

%     % xlabel('$\langle S/N \rangle$', 'interpreter', 'latex')
%     % ylabel('PDF')
%     % set(gca, 'FontSize',15)

%     % exportgraphics(fig, 'EW/hist_sigmaCIV_mSN.png', 'Resolution', 800)
errDR = abs(DR_1548_1550_DR12_flux).*sqrt((ErrREW_1548_flux./REW_1548_DR12_flux).^2 + ...
                                        (ErrREW_1550_flux./REW_1550_DR12_flux).^2);



end 

x = reshape(DR_1548_1550_DR12_flux(:,2:end), [], 1);
y = reshape(errDR(:,2:end),[],1);
pAll = reshape(p_c4(:,2:end), [],1);
nnz((x<1-y | x>2+y) & pAll>0.80)/nnz(pAll>0.80)
figure()

histogram(x(pAll>0.95))
xlim([0,3])