filename = 'data/dr7/processed/civ_samples_RJ-0-nAVG-1-num-Civ-50000-N-12.88-15.80-sgm-1-70-vCut-3000-alpha-90-fN-1.0e-01-fS-1.0e-01-Sgl-1-Multi-4.mat'; 
variables_to_load = {'offset_z_samples',...
                        'offset_sigma_samples', 'log_nciv_samples', 'nciv_samples'};

load(filename, variables_to_load{:});
fig =figure();
h = histogram(log_nciv_samples)
h.FaceColor = [0.3010 0.7450 0.9330];
h.NumBins=20;
set(get(gca, 'XLabel'), 'String', 'N_{CIV}');
set(get(gca, 'YLabel'), 'String', 'Frequency');
exportgraphics(fig, 'dist-N.png', 'Resolution', 400)