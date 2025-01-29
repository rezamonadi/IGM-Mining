function y=pltQSO(this_flux, this_wavelengths, c4_muL2, c4_muL1, ind_zoomL2, ind_zoomL1, z_EWlow, z_EWhigh, z_PM_test_plot,...
                   ind_not_remove, ttl, fid)
% function y =pltQSO(civ_1548_wavelength, this_flux, this_wavelengths, c4_muL2, c4_muL1, ...
%     this_sample_log_likelihoods_c4L2, sample_z_c4, log_nciv_samples, ttl,  fid)

    % been applied this_pixel_mask.
    
    fig = figure('visible', 'off');
    clf();
    % subplot('position', [0.05 0.49 0.90 5]);
    % construct dla_mu_map
    mgii_2796_wavelength= 2.7964e+03; 
    this_z_c4 = (this_wavelengths / mgii_2796_wavelength) - 1;

    % % subplot(2,1,1);
    
    % p = stairs(this_z_c4(ind_not_remove), this_flux(ind_not_remove));
    p = stairs(this_z_c4, this_flux);
    p.LineWidth = .5;
    p.Color = [0.3010 0.7450 0.9330, 0.8];
    hold on

    p = plot(this_z_c4, c4_muL2);
    p.Color = [0.8500 0.3250 0.0980, 0.6]; 
    p.LineWidth=1.5;
    hold on

    p = plot(this_z_c4, c4_muL1);
    p.Color = [0.500 0.8250 0.0980, 0.6]; 
    p.LineWidth=1.5;
    hold on 
    
    legend({'Flux', 'M$_D$', 'M$_S$'}, 'interpreter', 'latex')
    % legend({'Flux', 'M$_D$'}, 'interpreter', 'latex')
    hold on
    % xline(z_EWlow)
    % hold on
    % xline(z_EWhigh)
    % legend({'Flux', 'CIV', 'Singlet'}, 'interpreter', 'latex')
    xlim([min(this_z_c4), max(this_z_c4)])
    xlabel('$\lambda$/2796 (\AA) - 1', 'Interpreter','latex');
    ylabel('Normalised Flux');
    % title(ttl, 'FontSize', 5, 'interpreter', 'latex')

    for i=1:length(z_PM_test_plot)
        p=xline(z_PM_test_plot(i));
        p.Color = [0.1,0.1,0.1];
        p.LineStyle = '--';
        p.LineWidth=1;
        p.HandleVisibility = 'off';
        hold on 
    end
% %     
% %   
%     subplot('position', [0.05 0.49 0.90 5]);

%     subplot(2,1,2);
%     hold on
%     norm_sample_log_likelihoods = this_sample_log_likelihoods_c4L2;
%     norm_sample_log_likelihoods = norm_sample_log_likelihoods - nanmax(this_sample_log_likelihoods_c4L2);
%     norm_sample_log_likelihoods = norm_sample_log_likelihoods - log(sum(exp(norm_sample_log_likelihoods)));

    
%     s=scatter(sample_z_c4, log_nciv_samples, 5,...
%                     norm_sample_log_likelihoods, 'filled', 'DisplayName', 'sample likelihoods');
%     s.MarkerFaceAlpha = 0.5;
%     hcb = colorbar('southoutside');
%     caxis([0.5*nanmin(norm_sample_log_likelihoods) nanmax(norm_sample_log_likelihoods)]);
%     hcb.Label.String= 'log(Likelihood)';
%     xlabel('$z_{CIV}$', 'Interpreter','latex');
%     ylabel('$N_{CIV} (cm^{-2})$',  'Interpreter','latex');
%     xlim([min(this_z_c4), max(this_z_c4)])

    % axes('position', [0.2 0.75 0.18 0.15]);
    % box on 
    % p = plot(this_z_c4(ind_zoomL2), c4_muL2(ind_zoomL2));
    % p.Color = [0.8500 0.3250 0.0980, 0.6]; 
    % p.LineWidth=1;
    % hold on 
    % p = stairs(this_z_c4(ind_zoomL2), this_flux(ind_zoomL2));
    % p.Color = [0.3010 0.7450 0.9330];
    % hold on
    % xline(z_EWlow)
    % hold on
    % xline(z_EWhigh)
    title(ttl, 'FontSize', 7)

    % axes('position', [0.45 0.75 0.18 0.15]);
    % box on 
    % p = plot(this_z_c4(ind_zoomL1), c4_muL1(ind_zoomL1));
    % p.Color = [0.500 0.8250 0.0980, 0.6]; 
    % p.LineWidth=1;
    % hold on 
    % p = stairs(this_z_c4(ind_zoomL1), this_flux(ind_zoomL1));
    % p.Color = [0.3010 0.7450 0.9330];

    % set(gca, 'FontSize', 15);
    exportgraphics(fig, fid,'Resolution', 800)
    
end







