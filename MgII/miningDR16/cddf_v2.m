%---------- C IV CDDF ---------------
clear; close all; clc;
p_threshold =0.99;

% 1. LOAD DATA

% From posterior sampling
load('sampled_logN_per_absorber_P99_SNR4.mat');          

% From savingCat 
load('processed_qsos_dr16_parameters.mat');   

%catDR16 = load("catalog.mat");

% defining dX_dz function 

Omega_M = 0.3; Omega_Lambda = 0.7;
H0 = 70; c = 2.99792458e5;                     % km/s
H0_cgs = H0 * 1e5 / 3.085677581e24;
c_cgs = c * 1e5;                               % cm/s

dX_dz = @(z) (1+z).^2 ./ sqrt(Omega_M*(1+z).^3 + Omega_Lambda);

used_qsos = unique([sampledAbsorbers.qso_idx]);
z_min_all = savingCat.all_min_z_c4s(used_qsos);
z_max_all = savingCat.all_max_z_c4s(used_qsos);

% finding z_min and z_max for pathlines 
lambda_min_obs = 3600; lambda_max_obs = 10400; % SDSS spectral range
lambda_rest_c4 =  1548.1949462890625;	
z_c4_min_allowed = lambda_min_obs / lambda_rest_c4 - 1;
z_c4_max_allowed = lambda_max_obs / lambda_rest_c4 - 1;

z_min_clipped = max(z_min_all, z_c4_min_allowed);
z_max_clipped = min(z_max_all, z_c4_max_allowed);


%% Obtaining f(N,z)



ind_good = (SNR>SNR_cut);
nnz(ind_good)
z_abs_good = savingCat.all_map_z_c4L2;
z_abs_good = z_abs_good(ind_good);

log10_N_good = savingCat.all_map_N_c4L2;
log10_N_good = log10_N_good(ind_good);

h = histogram2(log10_N_good, z_abs_good, 'DisplayStyle','tile','ShowEmptyBins','on');
h.BinWidth = [0.2 0.2];

f_N_z = h.Values; 
Delta_log10_N = h.BinWidth(1);
log10_N_centers = h.XBinEdges(1:end-1) + h.BinWidth(1)*0.5;

% Telluric intervals in z for the transition lambda_rest_c4:
telluric_z = [ (6302-5)/lambda_rest_c4 - 1, (6302+5)/lambda_rest_c4 - 1;  % [zmin, zmax]
               (5579-5)/lambda_rest_c4 - 1, (5579+5)/lambda_rest_c4 - 1 ];

% helper: integrate dX/dz between a and b (0 if b<=a)
int_dX = @(a,b) (b>a) .* trapz(linspace(a,b,200), dX_dz(linspace(a,b,200)));

% --- ΔX_eff per z-bin: sum over sightlines, subtract tellurics per sightline ---
nSL = numel(z_min_clipped); % number of sightlines 
z_bin_centers = h.YBinEdges(1:end-1) + h.BinWidth(2)*0.5;
DeltaX_eff_per_zbin = zeros(1, numel(z_bin_centers));
for k = 1:numel(z_bin_centers)
    zb_lo = z_bin_centers(k) - 0.5*h.BinWidth(2);
    zb_hi = z_bin_centers(k) + 0.5*h.BinWidth(2);

    dX_sum = 0;
    for s = 1:nSL
        % overlap of sightline s with this z bin
        zlo = max(z_min_clipped(s), zb_lo);
        zhi = min(z_max_clipped(s), zb_hi);
        if zhi <= zlo, continue; end

        dX_s = int_dX(zlo, zhi);

        % subtract telluric overlaps (PER sightline)
        for t = 1:size(telluric_z,1)
            zt_lo = max(zlo, telluric_z(t,1));
            zt_hi = min(zhi, telluric_z(t,2));
            if zt_hi > zt_lo
                dX_s = dX_s - int_dX(zt_lo, zt_hi);
            end
        end

        dX_sum = dX_sum + dX_s;
    end

    DeltaX_eff_per_zbin(k) = dX_sum;
end


% --- Build the raw 2D counts (log10 N vs z) ---
% (Rebuild to ensure counts match h.BinWidth and edges)
h = histogram2(log10_N_good, z_abs_good, 'DisplayStyle','tile','ShowEmptyBins','on');
h.BinWidth = [0.2 0.2];
counts     = h.Values;
logN_ctr   = h.XBinEdges(1:end-1) + 0.5*h.BinWidth(1);
z_ctr      = h.YBinEdges(1:end-1) + 0.5*h.BinWidth(2);
Delta_logN = h.BinWidth(1);

% --- Normalize to φ(logN, z) = counts / (ΔlogN * ΔX_eff) ---
phi_logN_z = zeros(size(counts));
for k = 1:numel(z_ctr)
    dX_eff = DeltaX_eff_per_zbin(k);
    if dX_eff > 0 && isfinite(dX_eff)
        phi_logN_z(:,k) = counts(:,k) ./ (Delta_logN * dX_eff);
    else
        phi_logN_z(:,k) = 0;
    end
end


% --- Plot φ(logN, z) as heatmap ---
fig = figure;
histogram2('XBinEdges', h.XBinEdges, 'YBinEdges', h.YBinEdges, 'BinCounts', phi_logN_z, ...
           'DisplayStyle','tile','ShowEmptyBins','on');
colorbar; xlabel('log_{10} N'); ylabel('z'); title('$\phi(\log N, z)$  [per dlogN per dX]', 'Interpreter','latex');
exportgraphics(fig, 'figsSN4/Phi-2D-P99_SNR4.png', 'Resolution', 1000)

% --- Redshift-averaged 1D CDDF over your full z-range (path-length weighted) ---
X_tot  = sum(DeltaX_eff_per_zbin);
phi_1D = (phi_logN_z * DeltaX_eff_per_zbin.') / X_tot;  % size = nLogN x 1


fig = figure;
plot(logN_ctr, phi_1D./(log(10)*(10.^logN_ctr.')), '-o')
xlabel('$log_{10}$ N', 'Interpreter', 'latex')
ylabel('CDDF')
grid on
title('Redshift-averaged CDDF')
yscale('log')
exportgraphics(fig, 'figsSN4/CDDF-Avg-P99_SNR4.png', 'Resolution', 1000)

% line density 
fig=figure;

dNdX_per_z = sum(phi_logN_z, 1) * Delta_logN;   % [1 x nZ]
plot(z_ctr, (dNdX_per_z)); yscale('log')
exportgraphics(fig, 'figsSN4/lineDensity-P99_SNR4.png',  'Resolution', 1000)


