ew = load('REW_1548_DR12.mat');
ewv = ew.REW_1548_DR12_voigt;
ewf = ew.REW_1548_DR12_flux;

fig = figure();           
DREW = ewv - ewf;
DREW = reshape(DREW, length(ewv)*7, 1);
DREW_cut = DREW(DREW>quantile(DREW, 0.01) & DREW<quantile(DREW, 0.99));
p = histogram(DREW_cut);
p.NumBins=50;
xlabel('REW(CIV)\_Voigt - REW(CIV)\_Flux')
ylabel('Abundance')
exportgraphics(fig, 'DREWDR12.png', 'Resolution', 800);

fig = figure();           
rREW = 1 - ewv./ewf;
rREW = reshape(rREW, length(ewv)*7, 1);
% rREW_cut = DREW(rREW>quantile(rREW, 0.01) & rREW<quantile(rREW, 0.99));
rREW_cut = rREW;
p = histogram(rREW_cut);
p.NumBins=50;
xlabel('1 - REW(CIV)\_Voigt/REW(CIV)\_Flux')
ylabel('Abundance')
exportgraphics(fig, 'rREWDR12.png', 'Resolution', 800);