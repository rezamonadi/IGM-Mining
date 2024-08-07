Ndr12 = reshape(all_map_N_c4L2(all_p_c4>0.95), [],1);
z = reshape(all_map_z_c4L2(all_p_c4>0.95), [],1);


Ndr7 = all_N_civ;
Ndr7 = Ndr7(Ndr7>0);
Ndr7 = reshape(Ndr7,[],1);
z_dr7 = reshape(all_z_civ3, [], 1);
z_dr7 = z_dr7(Ndr7>0);
fig = figure();
t = tiledlayout(2, 3);

ax1 = nexttile;
hh = histogram(ax1, Ndr12(z<=1.5), 'numbins',20, 'normalization', 'PDF');
binEdges = hh.BinEdges;
hold on 
hh2 = histogram(ax1, Ndr7(z_dr7<=1.5), 'normalization', 'pdf');
hh2.BinEdges = binEdges;
hold on
% xline(median(Ndr12(z<=1.5)), 'linestyle', '--')
hold on
text(ax1, 15, 1.2, 'z<1.5')


ax2 = nexttile;
hh = histogram(ax2, Ndr12(z>1.5 & z<=2), 'numbins',20, 'normalization', 'pdf');
binEdges = hh.BinEdges;
hold on
hh2 = histogram(ax2, Ndr7(z_dr7>1.5 & z_dr7 <=2), 'normalization', 'PDF');
hh2.BinEdges = binEdges;
% xline(median(Ndr12(z>1.5 & z<=2)), 'linestyle', '--')
text(ax2, 14.7, 1.2, '1.5>z>2')


ax3 = nexttile;
hh = histogram(ax3, Ndr12(z>2 & z<=2.5), 'numbins',20, 'normalization', 'pdf');
binEdges = hh.BinEdges;
hold on
hh2 = histogram(ax3, Ndr7(z_dr7>2 & z_dr7<=2.5), 'normalization','pdf');
hh2.BinEdges = binEdges;
% xline(median(Ndr12(z>2 & z<=2.5)), 'linestyle', '--')
text(ax3, 14.7, 1.2, '2>z>2.5')

ax4 = nexttile;
hh  = histogram(ax4, Ndr12(z>2.5 & z<=3), 'numbins',20, 'normalization', 'pdf');
binEdges = hh.BinEdges;
hold on
hh2 = histogram(ax4, Ndr7(z_dr7>2.5 & z_dr7<=3), 'normalization','pdf');
hh2.BinEdges = binEdges;
text(ax4, 14.8, 1.2, '2.5>z>3')
% xline(median(Ndr12(z>2.5 & z<=3)), 'linestyle', '--')

ax5 = nexttile;
hh = histogram(ax5, Ndr12(z>3 & z<=3.5), 'numbins',20, 'normalization', 'pdf');
binEdges = hh.BinEdges;
hold on
hh2 = histogram(ax5, Ndr7(z_dr7>3 & z_dr7<=3.5), 'normalization','pdf');
hh2.BinEdges = binEdges;
% xline(median(Ndr12(z>3 & z<=3.5)), 'linestyle', '--')
text(ax5, 14.8, 1.2, '3>z>3.5')

ax6 = nexttile;
hh = histogram(ax6, Ndr12(z>3.5), 'numbins',20, 'normalization', 'pdf')
binEdges = hh.BinEdges;
hold on
% xline(median(Ndr12(z>3.5)), 'linestyle', '--')
text(ax6, 15, .95, 'z>3.5')
hh2 = histogram(ax6, Ndr7(z_dr7>3.5), 'normalization','pdf')
hh2.BinEdges = binEdges;
t.Padding = 'compact';
t.TileSpacing = 'compact';

xlabel(t, 'log(N_{CIV})');
ylabel(t, 'PDF')

exportgraphics(fig, 'NDR12-zBands.png', 'resolution',800)