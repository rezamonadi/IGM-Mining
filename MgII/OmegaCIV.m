clear
load('ColumnDensityDistRedShiftBinsData.mat')
numRedShiftBins = 8;
BinsRedshift = linspace(1, 5, numRedShiftBins + 1);
widthRedShiftBins = 0.5;
midBinsRedShift = BinsRedshift(2:end) - 0.5*widthRedShiftBins;

for j=1:18
	AllMeanDensityBin(j) = mean(AllCountColumnDensity(:,j)); 
	for k=1:numRedShiftBins
		AllMeanDensityRedShift(j, k) = mean(AllCountColumnDensityRedShift(:, j, k));
	end
end

fig = figure();
semilogy(midBinColumnDensity', AllMeanDensityBin(:), 'lineWidth', 3, 'Color', [0.1 0.1 0.1]) 
xlabel('$log({\rm N}_{\rm CIV})$', 'interpreter', 'latex')
ylabel('Count')
exportgraphics(fig, 'N_unc.png', 'Resolution', 800)


fig = figure();
for k=1:numRedShiftBins
	semilogy(midBinColumnDensity', AllMeanDensityRedShift(:,k))
	hold on 
	
end
legend('1<z<1.5', '1.5<z=2', '2<z<2.5', '2.5<z<3', '3<z<3.5', '3.5<z<4', '4<z<4.5', '4.5<z<5')
exportgraphics(fig, 'N_zBin.png', 'Resolution', 800)
		

Omega_M = 0.285;
Omega_Lambda = 0.742;
dX = @(z) (1+z).^2./(sqrt(Omega_M*(1+z).^3 + Omega_Lambda)/(3*Omega_M)); 
Delta_X = @(z) (integral(dX, 0, z));
numRedShiftBins = 8;
BinsRedshift = linspace(1, 5, numRedShiftBins + 1);
widthRedShiftBins = 0.5;
midBinsRedShift = BinsRedshift(2:end) - 0.5*widthRedShiftBins;

ProtonMass = 2e-23; %proton mass in g
c   = 2.99792458e+10;   % cm/s
h100 = 3.2407789e-18 * 0.719;
rho_c = 9.77e-30;
conv = ProtonMass*h100/c/rho_c;
for i=1:numRedShiftBins
    Omega_civ(i) = conv*sum(10.^midBinColumnDensity'.*AllMeanDensityRedShift(i))./(Delta_X(BinsRedshift(i+1)) - Delta_X(BinsRedshift(i)));
    Omega_civ14p(i) = conv*sum(10.^midBinColumnDensity(8:end)'.*AllMeanDensityRedShift(8:end, i))./(Delta_X(BinsRedshift(i+1)) - Delta_X(BinsRedshift(i)));
end

fig = figure();
plot(midBinsRedShift, log(Omega_civ*1e8))
hold on
plot(midBinsRedShift, log(Omega_civ14p*1e8))
ylabel('$\Omega_{\rm CIV} (10^{-8})$', 'interpreter', 'latex')
xlabel('$Z_{\rm CIV}$', 'interpreter', 'latex')
exportgraphics(fig, 'Omega_civN14Plus.png', 'Resolution', 800)

