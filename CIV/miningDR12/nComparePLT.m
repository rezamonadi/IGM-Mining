close all
clear 
load('ShortProcessedDR12.mat');
set_parameters_dr7;
build_catalog_dr7;

zCIV_PM = all_z_civ3;              
nCIV_PM = all_N_civ;
nCIV_GP = all_map_N_c4L2;
zCIV_GP = all_map_z_c4L2;
pCIV = all_p_c4;

nCIV_GP_1D = reshape(nCIV_GP, [], 1);
zCIV_GP_1D = reshape(zCIV_GP, [],1);
pCIV_1D    = reshape(pCIV, [], 1);

nCIV_GP_1D_01 = nCIV_GP_1D(pCIV_1D>0.01);
zCIV_GP_1D_01 = zCIV_GP_1D(pCIV_1D>0.01);

nCIV_GP_1D_65 = nCIV_GP_1D(pCIV_1D>0.65);
zCIV_GP_1D_65 = zCIV_GP_1D(pCIV_1D>0.65);

nCIV_GP_1D_95 = nCIV_GP_1D(pCIV_1D>0.95);
zCIV_GP_1D_95 = zCIV_GP_1D(pCIV_1D>0.95);

nCIV_GP_1D_99 = nCIV_GP_1D(pCIV_1D>0.99);
zCIV_GP_1D_99 = zCIV_GP_1D(pCIV_1D>0.99);




nCIV_PM_1D = reshape(nCIV_PM, [], 1);
zCIV_PM_1D = reshape(zCIV_PM, [],1);

nCIV_PM_1D_c = nCIV_PM_1D(zCIV_PM_1D>0 & nCIV_PM_1D>0);
zCIV_PM_1D_c = zCIV_PM_1D(zCIV_PM_1D>0 & nCIV_PM_1D>0);



h1 = histogram(nCIV_GP_1D_01);
h1.BinWidth=0.2;
bins = h1.BinEdges;
x_plot = bins(2:end) - h1.BinWidth;
y1 = h1.Values;

h2 = histogram(nCIV_GP_1D_65);
h2.BinEdges = bins;
y2 = (h2.Values)./y1;

h3 = histogram(nCIV_GP_1D_95);
h3.BinEdges = bins;
y3 = (h3.Values)./y1;

h4 = histogram(nCIV_GP_1D_99);
h4.BinEdges = bins;
y4 = (h4.Values)./y1;

close all
fig = figure();

plot(x_plot, y2, 'LineWidth', 3)
hold on
plot(x_plot, y3, 'LineWidth', 3)
hold on
plot(x_plot, y4, 'LineWidth', 3)
ylabel('Count/Count(P(M_D)>0.01)')
xlabel('$\rm{N}_{\rm CIV}$', 'interpreter','latex')
legend('P(M_D)>0.65', 'P(M_D)>0.95', 'P(M_D)>0.99')
exportgraphics(fig, 'N_compare.png', 'Resolution', 800)










% fig = figure()
% scatter(zCIV_GP_1D, nCIV_GP_1D, 10, 'filled', 'MarkerFaceAlpha', 0.5)
% hold on 
% scatter(zCIV_PM_1D_c, nCIV_PM_1D_c, 10, 'filled', 'MarkerFaceAlpha', 0.5)
% xlabel('$\rm{z}_{\rm CIV}$', 'interpreter','latex')
% ylabel('$\rm{N}_{\rm CIV}$', 'interpreter','latex')
% legend({'GP(P(D)>0.95)','PM'})
% exportgraphics(fig, 'nCompare_PM_GP.png', 'resolution',400)


% fig = figure()
% XY_bin = bin_averager(zCIV_GP_1D,nCIV_GP_1D, 7);
% plot(XY_bin(:,1), XY_bin(:,2))
   
% hold on 
% XY_bin = bin_averager(zCIV_PM_1D_c,nCIV_PM_1D_c, 7);
% plot(XY_bin(:,1), XY_bin(:,2))
% xlabel('$\rm{z}_{\rm CIV}$', 'interpreter','latex')
% ylabel('$\rm{N}_{\rm CIV}$', 'interpreter','latex')
% legend({'GP(P(D)>0.95)','PM'}, 'Location', 'northwest')
% exportgraphics(fig, 'nCompare_PM_GP_line_avg_bin7.png', 'resolution',400)






