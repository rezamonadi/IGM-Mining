function y=plt_mu(this_flux, this_wavelengths, mu, z_qso, ttl, fid)


    % been applied this_pixel_mask.
    
    fig = figure('visible', 'off');
    clf();
    % subplot('position', [0.05 0.49 0.90 1.45]);
    % construct dla_mu_map
    % this_z_c4 = (this_wavelengths / 1549.48) - 1;

    % subplot(2,1,1);
    this_rest_wavelengths = this_wavelengths/(1+z_qso);
    p = stairs(this_rest_wavelengths, this_flux);
    p.LineWidth = .3;
    p.Color = [0.1 0.2 0.9, 0.8];
    hold on

    

    p = plot(this_rest_wavelengths, mu);
    p.Color = [0.900 0.1 0.1, 0.5]; 
    p.LineWidth=2;
    
    % legend( 'Flux', 'Continuumn')
    hold on
    xlim([min(this_rest_wavelengths), max(this_rest_wavelengths)]+10)
    % xlabel('(Observed Wavelengths $\lambda$ (\AA) / 1549.48 (\AA)) - 1', 'FontSize', 16, 'Interpreter','latex');
    xlabel('$\lambda_r$ (\AA)', 'FontSize', 16, 'Interpreter','latex');
    ylabel('Normalized Flux ',                            'FontSize', 16, 'Interpreter','latex');
    title(ttl, 'FontSize', 14)
    
    CIV = 1548.2040;
    CII = 1335.31;
    SIV = 1397.61;
    vl = xline(CIV, '--',  'CIV'); 
    vl.FontSize = 20;
    vl = xline(SIV, '--',  'SIV'); 
    vl.FontSize = 20;
    vl = xline(CII, '--',  'CII'); 
    vl.FontSize = 20;
    
    set(gca,'FontSize',20);
    


    exportgraphics(fig, fid,'Resolution', 800)
    
end







