function y = pltQSOv1(this_flux, this_wavelengths, c4_muL2, ttl, fid)
%PLTQSOV1 Plot normalised flux/model and their residual on aligned panels.
%   The residual is Flux - M_D, where M_D is the model drawn in the top panel.

    civ_1548_wavelength = 1548.1949462890625;
    this_z_c4 = (this_wavelengths / civ_1548_wavelength) - 1;
    % residual = this_flux - c4_muL2;

    fig = figure('visible', 'off');
    tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax_flux = nexttile;
    p_flux = stairs(ax_flux, this_z_c4, this_flux);
    p_flux.LineWidth = 0.5;
    p_flux.Color = [0.3010 0.7450 0.9330 0.8];
    hold(ax_flux, 'on');

    p_model = plot(ax_flux, this_z_c4, c4_muL2);
    p_model.Color = [0.8500 0.3250 0.0980 0.6];
    p_model.LineWidth = 1.5;
    legend(ax_flux, {'Flux', 'M$_D$'}, 'Interpreter', 'latex');
    ylabel(ax_flux, 'Normalised Flux');
    title(ax_flux, ttl, 'FontSize', 7);
    xlim(ax_flux, [min(this_z_c4), max(this_z_c4)]);
    ax_flux.XTickLabel = [];

    ax_residual = nexttile;
    p_residual = stairs(ax_residual, this_z_c4, residual);
    p_residual.LineWidth = 0.5;
    p_residual.Color = [0.3010 0.7450 0.9330 0.8];
    hold(ax_residual, 'on');
    yline(ax_residual, 0, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.75, ...
          'HandleVisibility', 'off');
    xlim(ax_residual, [min(this_z_c4), max(this_z_c4)]);
    xlabel(ax_residual, '$\lambda$/1548 (\AA) - 1', 'Interpreter', 'latex');
    ylabel(ax_residual, 'Flux - M$_D$', 'Interpreter', 'latex');
    linkaxes([ax_flux, ax_residual], 'x');

    exportgraphics(fig, fid, 'Resolution', 800);
    y = [];
end
