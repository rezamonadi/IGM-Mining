function y=pltQSO(this_flux, this_wavelengths, c4_muL2, c4_muL1, ind_zoomL2, ind_zoomL1, ttl, fid)


    % been applied this_pixel_mask.
    
    fig = figure('visible', 'off');
    clf();
    % subplot('position', [0.05 0.49 0.90 1.45]);
    % construct dla_mu_map
    this_z_c4 = (this_wavelengths / 1550) - 1;

    % subplot(2,1,1);
    
    p = stairs(this_z_c4, this_flux);
    p.LineWidth = .2;
    p.Color = [0.3010 0.7450 0.9330];
    
    hold on
    p = plot(this_z_c4, c4_muL2);
    p.Color = [0.8500 0.3250 0.0980, 0.6]; 
    p.LineWidth=1;
    hold on

    % hold on
    p = plot(this_z_c4, c4_muL1);
    p.Color = [0.500 0.8250 0.0980, 0.6]; 
    p.LineWidth=1;

    legend( 'Flux', 'C4', 'S')
    % legend( 'Flux', 'C4')
    xlim([min(this_z_c4), max(this_z_c4)])
    xlabel('(observed wavelengths $\lambda$ (\AA) / 1549 (\AA)) - 1', 'FontSize', 14, 'Interpreter','latex');
    ylabel('normalized flux $\mathbf{y}$',                            'FontSize', 14, 'Interpreter','latex');
    title(ttl, 'FontSize', 7)
    
   
    % subplot(2,1,2);
    % hold on
    % norm_sample_log_likelihoods = this_sample_log_likelihoods_c4L2;
    % norm_sample_log_likelihoods = norm_sample_log_likelihoods - nanmax(this_sample_log_likelihoods_c4L2);
    % norm_sample_log_likelihoods = norm_sample_log_likelihoods - log(sum(exp(norm_sample_log_likelihoods)));

    
    % s=scatter(sample_z_c4, log_nciv_samples, 5,...
    %                 norm_sample_log_likelihoods, 'filled', 'DisplayName', 'sample likelihoods');
    % s.MarkerFaceAlpha = 0.5;
    % hcb = colorbar('southoutside');
    % caxis([0.5*nanmin(norm_sample_log_likelihoods) nanmax(norm_sample_log_likelihoods)]);
    % hcb.Label.String= 'log(Likelihood)';
    % xlabel('$z_{CIV}$', 'FontSize', 10, 'Interpreter','latex');
    % ylabel('$\log N_{CIV}$', 'FontSize', 10, 'Interpreter','latex');
    xlim([min(this_z_c4), max(this_z_c4)])

    axes('position', [0.2 0.68 0.18 0.15]);
    box on 
    p = plot(this_z_c4(ind_zoomL2), c4_muL2(ind_zoomL2));
    p.Color = [0.8500 0.3250 0.0980, 0.6]; 
    p.LineWidth=1;
    hold on 
    p = stairs(this_z_c4(ind_zoomL2), this_flux(ind_zoomL2));
    p.Color = [0.3010 0.7450 0.9330];

    axes('position', [0.45 0.68 0.18 0.15]);
    box on 
    p = plot(this_z_c4(ind_zoomL1), c4_muL1(ind_zoomL1));
    p.Color = [0.500 0.8250 0.0980, 0.6]; 
    p.LineWidth=1;
    hold on 
    p = stairs(this_z_c4(ind_zoomL1), this_flux(ind_zoomL1));
    p.Color = [0.3010 0.7450 0.9330];


    exportgraphics(fig, fid,'Resolution', 800)
    return 
end







