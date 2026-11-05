function y = pltQSOv2(this_flux, this_wavelengths, c4_muL2, c4_muL1, ttl, fid)
%PLTQSOV2 Plot flux with the current doublet and singlet models.
%   c4_muL2 and c4_muL1 are model flux vectors evaluated on the same
%   surviving pixels as this_flux and this_wavelengths.

    civ_1548_wavelength = 1548.1949462890625;
    this_z_c4 = (this_wavelengths / civ_1548_wavelength) - 1;

    fig = figure('visible', 'off');
    figure_cleanup = onCleanup(@() close(fig));
    tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax_doublet = nexttile;
    p_flux = stairs(ax_doublet, this_z_c4, this_flux);
    p_flux.LineWidth = 0.5;
    p_flux.Color = [0.3010 0.7450 0.9330 0.8];
    hold(ax_doublet, 'on');
    p_model = plot(ax_doublet, this_z_c4, c4_muL2);
    p_model.Color = [0.8500 0.3250 0.0980 0.6];
    p_model.LineWidth = 1.5;
    legend(ax_doublet, {'Flux', 'M$_D$'}, 'Interpreter', 'latex');
    ylabel(ax_doublet, 'Normalised Flux');
    title(ax_doublet, ttl, 'FontSize', 7);
    xlim(ax_doublet, [min(this_z_c4), max(this_z_c4)]);
    ax_doublet.XTickLabel = [];

    ax_singlet = nexttile;
    p_flux = stairs(ax_singlet, this_z_c4, this_flux);
    p_flux.LineWidth = 0.5;
    p_flux.Color = [0.3010 0.7450 0.9330 0.8];
    hold(ax_singlet, 'on');
    p_model = plot(ax_singlet, this_z_c4, c4_muL1);
    p_model.Color = [0.8500 0.3250 0.0980 0.6];
    p_model.LineWidth = 1.5;
    legend(ax_singlet, {'Flux', 'M$_S$'}, 'Interpreter', 'latex');
    ylabel(ax_singlet, 'Normalised Flux');
    xlabel(ax_singlet, '$\lambda$/1548 (\AA) - 1', 'Interpreter', 'latex');
    xlim(ax_singlet, [min(this_z_c4), max(this_z_c4)]);
    common_ylim = [min([ax_doublet.YLim, ax_singlet.YLim]), ...
                   max([ax_doublet.YLim, ax_singlet.YLim])];
    ylim(ax_doublet, common_ylim);
    ylim(ax_singlet, common_ylim);
    linkaxes([ax_doublet, ax_singlet], 'xy');

    exportgraphics(fig, fid, 'Resolution', 800);
    y = [];
end
