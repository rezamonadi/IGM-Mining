function y=plt_mu(this_flux, this_wavelengths, mu, z_qso, M, ttl, fid)

    set_parameters_dr7
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
    p.LineWidth=1;
        hold on 


    DD = diag(M*M');
    p = plot(this_rest_wavelengths, mu + DD);
    p.Color = [0.900 0.1 0.1, 0.5]; 
    p.LineWidth = 0.5;
    hold on 


    p = plot(this_rest_wavelengths, mu - DD);
    

    p.Color = [0.900 0.1 0.1, 0.5]; 
    p.LineWidth = 0.5;
    % hold on
    X = [this_rest_wavelengths', fliplr(this_rest_wavelengths')];
    Y = [mu'+DD', fliplr(mu'-DD')];
    p = fill(X,Y, 'r');
    p.FaceAlpha = 0.3;
    % legend( 'Flux', 'Continuum')
    hold on
   
    % vl = xline(max(this_rest_wavelengths) - kms_to_z(vCut)*mgii_2796_wavelength);
    % vl.LineWidth = 1
    % vl.Color = [1, 0.1, 0.2];
    % hold on
    % y_points = [min(this_flux)-0.2,max(this_flux), max(this_flux), min(this_flux)-0.2]
    % color = [1,1,0]
    % a = fill(x_points, y_points, color)
    % a.FaceAlpha = 0.1;
    % hold on
    xlim([min(this_rest_wavelengths), max(this_rest_wavelengths)]+10)
    xlabel('Rest Wavelengths (\AA)', 'FontSize', 16, 'Interpreter','latex');
    % xlabel('$\lambda_r$ (\AA)', 'FontSize', 16, 'Interpreter','latex');
    ylabel('Normalized Flux ', 'FontSize', 16, 'Interpreter','latex');
    % title(ttl, 'FontSize', 14)
    ylim([min(this_flux)-0.1, max(mu+DD)+0.1])
    CIV = 1548.2040;
    CII = 1335.31;
    SIV = 1397.61;
    CIIIF = 1908.734;
    MgII = 2799.117  ;
    vl = xline(CIV, '--',  'CIV'); 
    vl.FontSize = 15;
    vl = xline(SIV, '--',  'SIV'); 
    vl.FontSize = 15;
    vl = xline(CII, '--',  'CII'); 
    vl.FontSize = 15;
    vl = xline(CIIIF, '--',  'CIII]'); 
    vl.FontSize = 15;
    vl = xline(MgII, '--',  'MgII'); 
    vl.FontSize = 15;
    set(gca,'FontSize',17);
    


    exportgraphics(fig, fid,'Resolution', 800)
    
end