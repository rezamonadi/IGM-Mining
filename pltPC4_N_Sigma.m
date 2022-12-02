fprintf('loading Processed File ...\n');
fname = 'data/dr12/processed/processed_dr12.mat';
load(fname)  

% loading full 3D poster
map_N_c4L2     = savingCat.all_map_N_c4L2;
map_sigma_c4L2 = savingCat.all_map_sigma_c4L2;  
p_c4           = savingCat.all_p_c4; 

test_ind       = savingCat.test_ind;
nqsos = nnz(test_ind);
AllSigma = reshape(map_sigma_c4L2, nqsos*7,1);
AllN = reshape(map_N_c4L2, nqsos*7,1);
AllPc4 = reshape(p_c4,  nqsos*7,1);
mask = AllPc4>0.5;
fig = figure();
p = scatter(AllN(mask), AllSigma(mask), 7, AllPc4(mask));
p.MarkerFaceColor = 'flat';
p.MarkerFaceAlpha = 0.3;
h = colorbar;
ylabel(h,'P(CIV)')
set(get(gca, 'XLabel'), 'String', 'N_{CIV}');
set(get(gca, 'YLabel'), 'String', '\sigma_{CIV}');
xlim([13.25 15])
exportgraphics(fig, 'P50NSigma.png', 'Resolution', 800)


