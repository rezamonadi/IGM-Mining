%---------- C IV CDDF ---------------
clearvars; close all; clc;

%% loading data
load('sampled_logN_per_absorber_Smp1000.mat');   % contains sampledAbsorbers struct array
load('processed_qsos_dr16_parameters.mat');      % contains savingCat struct array
catDR16 = load('catalog.mat');

%% 
SNR_threshold = 4;
p_threshold = 0.95; 
% Optional mask (uncomment p_c4 threshold if you want it applied here)
mask = ([sampledAbsorbers.SNR] > SNR_threshold) & ([sampledAbsorbers.p_c4] > p_threshold);
sampledAbsorbers = sampledAbsorbers(mask);

% Vectorized field extraction (much faster than loops)
p_c4 = [sampledAbsorbers.p_c4]';                           % n×1

% logN_samples: assume each element has a 1×100 row
logN_samples = vertcat(sampledAbsorbers.logN_samples);     % n×100

% z_abs_samples:
% If it's a scalar per absorber, replicate across 100 columns.
% If it's already 1×100 per absorber, just stack like logN_samples.
if isscalar(sampledAbsorbers(1).z_abs_samples)
    z_abs_samples = repmat([sampledAbsorbers.z_abs_samples]', 1, size(logN_samples,2));  % n×100
else
    z_abs_samples = vertcat(sampledAbsorbers.z_abs_samples);                              % n×100
end

% soft HistoGram --> Obtaining phi(log10(N),z) with soft using sampled logN values and fixed MAP(z_abs) values
log10_N_edges = linspace(min(min(logN_samples)), max(max(logN_samples)), 20);  % 20 bins from 12 to 16
z_edges       =  linspace(min(min(z_abs_samples)), max(max(z_abs_samples)), 20);
Out = soft_hist_2d(p_c4, logN_samples, z_abs_samples, log10_N_edges, z_edges);

%% % defining dX_dz function 

Omega_M = 0.3; Omega_Lambda = 0.7;
H0 = 70; c = 2.99792458e5;                     % km/s
H0_cgs = H0 * 1e5 / 3.085677581e24;
c_cgs = c * 1e5;                               % cm/s

dX_dz = @(z) (1+z).^2 ./ sqrt(Omega_M*(1+z).^3 + Omega_Lambda);


used_qsos = unique([sampledAbsorbers.qso_idx]);
z_min_all = savingCat.all_min_z_c4s(used_qsos);
z_max_all = savingCat.all_max_z_c4s(used_qsos);
% % finding z_min and z_max for pathlines 
lambda_min_obs = 3600; lambda_max_obs = 10400; % SDSS spectral range
lambda_rest_c4 =  1548.1949462890625;	
z_c4_min_allowed = lambda_min_obs / lambda_rest_c4 - 1;
z_c4_max_allowed = lambda_max_obs / lambda_rest_c4 - 1;

z_min_clipped = max(z_min_all, z_c4_min_allowed);
z_max_clipped = min(z_max_all, z_c4_max_allowed);


Delta_log10_N = Out.dx(1);  % assuming uniform bin width
log10_N_centers = Out.xcenters; % assuming uniform bin width

% Telluric intervals in z for the transition lambda_rest_c4:
telluric_z = [ (6302-5)/lambda_rest_c4 - 1, (6302+5)/lambda_rest_c4 - 1;  % [zmin, zmax]
                (5579-5)/lambda_rest_c4 - 1, (5579+5)/lambda_rest_c4 - 1 ];

% helper: integrate dX/dz between a and b (0 if b<=a)
int_dX = @(a,b) (b>a) .* trapz(linspace(a,b,200), dX_dz(linspace(a,b,200)));

% --- ΔX_eff per z-bin: sum over sightlines, subtract tellurics per sightline ---
nSL = numel(z_min_clipped); % number of sightlines 
z_bin_centers = Out.ycenters; % assuming uniform bin width
DeltaX_eff_per_zbin = zeros(1, numel(z_bin_centers));
for k = 1:numel(z_bin_centers)
    zb_lo = z_bin_centers(k) - 0.5*Out.dy(1);
    zb_hi = z_bin_centers(k) + 0.5*Out.dy(1);

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


%% --- Build the raw 2D counts (log10 N vs z) ---
% (Rebuild to ensure counts match h.BinWidth and edges)
counts     = Out.H2D;  % size = nLogN x nZ
err_counts = Out.sigma2D; 
logN_ctr   = Out.xcenters; % assuming uniform bin width
z_ctr      = Out.ycenters; % assuming uniform bin width
Delta_logN = Out.dx(1);  % assuming uniform bin width
Delta_z    = Out.dy(1);  % assuming uniform bin width

% --- Normalize to φ(logN, z) = counts / (ΔlogN * ΔX_eff) ---
phi_logN_z = zeros(size(counts));
err_phi_logN_z = zeros(size(counts));
for k = 1:numel(z_ctr)
    dX_eff = DeltaX_eff_per_zbin(k);
    if dX_eff > 0 && isfinite(dX_eff)
        phi_logN_z(:,k) = counts(:,k) ./ (Delta_logN * dX_eff);
        err_phi_logN_z(:,k) = err_counts(:,k) ./ (Delta_logN * dX_eff);

    else
        phi_logN_z(:,k) = 0;
        err_phi_logN_z(:,k) = 0;
    end
end


%% --- Plot φ(logN, z) as heatmap ---
fig = figure;
imagesc(logN_ctr, z_ctr, phi_logN_z.' ./ (log(10)*(10.^logN_ctr.'))); % convert to per dex
colorbar;
xlabel('$log_{10}(N_{CIV})$', 'Interpreter', 'latex');
ylabel('$z_{\rm{CIV}}$', 'Interpreter', 'latex');
title('$\phi(log_{10}(N_{CIV}), z)$  [per dlogN per dX]', 'Interpreter','latex');
exportgraphics(fig, 'figSoftP95snr4/phi.png', 'Resolution', 400)

%% --- Redshift-averaged 1D CDDF over your full z-range (path-length weighted) ---
X_tot  = sum(DeltaX_eff_per_zbin);
phi_1D     = (phi_logN_z.'     ./ (log(10)*(10.^logN_ctr.')) * DeltaX_eff_per_zbin.') / X_tot;  % size = nLogN x 1
err_phi_1D = (err_phi_logN_z.' ./ (log(10)*(10.^logN_ctr.')) * DeltaX_eff_per_zbin.') / X_tot;  % size = nLogN x 1
fig = figure;
errorbar(logN_ctr, phi_1D, err_phi_1D , '-o');
xlabel('$\log_{10}(N_{CIV})$', 'Interpreter', 'latex')
ylabel('$\frac{d\mathcal{N}}{d\mathcal{X}_{eff}}$', 'Interpreter','latex')
grid on
title('Redshift-averaged CDDF')
yscale('log')
exportgraphics(fig, 'figSoftP95snr4/phi1D.png', 'Resolution', 400)


%% line density 
fig=figure;

dNdX_per_z = sum(phi_logN_z, 1) * Delta_logN;   % [1 x nZ]
errorbar(z_ctr, (dNdX_per_z), Out.sigmaHy); yscale('log')
exportgraphics(fig, sprintf('figSoftP95snr4/lineDensity-P%2d.png', p_threshold*100), 'Resolution', 1000)

%% CDDF for custom z-ranges (exact ΔX weighting, no linear-in-z assumption)
nZ      = numel(z_edges) - 1;

% ---- redshift ranges (rows = [z_min z_max]) ----
z_ranges = [0.5 1; 1.0 2.0; 2.0 3.0; 3.0 5.0];

nRanges = size(z_ranges,1);
phi_1D_multi = zeros(numel(logN_ctr), nRanges);

% W_range(r,j): effective path length (ΔX) contributed by z-bin j inside z-range r
W_range = zeros(nRanges, nZ);

for r = 1:nRanges
    z_lo_r = z_ranges(r,1);
    z_hi_r = z_ranges(r,2);

    for j = 1:nZ
        % intersection of this z-bin with the requested z-range
        zb_lo = z_edges(j);
        zb_hi = z_edges(j+1);
        lo = max(zb_lo, z_lo_r);
        hi = min(zb_hi, z_hi_r);
        if hi <= lo, continue; end

        % sum ΔX over sightlines within [lo, hi], subtract tellurics per sightline
        dX_sum = 0;
        for s = 1:nSL
            % overlap with sightline s coverage
            zlo = max(z_min_clipped(s), lo);
            zhi = min(z_max_clipped(s), hi);
            if zhi <= zlo, continue; end

            dX_s = int_dX(zlo, zhi);

            % subtract telluric overlaps for this sightline segment
            for t = 1:size(telluric_z,1)
                zt_lo = max(zlo, telluric_z(t,1));
                zt_hi = min(zhi, telluric_z(t,2));
                if zt_hi > zt_lo
                    dX_s = dX_s - int_dX(zt_lo, zt_hi);
                end
            end

            dX_sum = dX_sum + dX_s;
        end

        W_range(r,j) = dX_sum;  % ΔX contributed by this z-bin inside the z-range
    end

    X_range = sum(W_range(r,:));         % total ΔX in this z-range
    if X_range > 0
        % Path-length–weighted average of φ(logN,z) over only the ΔX that lies in-range
        phi_1D_multi(:,r) = (phi_logN_z.' ./ (log(10)*(10.^logN_ctr.')) * W_range(r,:).') / X_range;   % [nLogN x 1]
    else
        phi_1D_multi(:,r) = 0;
        warning('No path in z-range [%.2f, %.2f]', z_lo_r, z_hi_r);
    end
end

% Plot
fig = figure; hold on; box on; grid on;
for r = 1:nRanges
    if all(phi_1D_multi(:,r)==0), continue; end
    plot(logN_ctr, phi_1D_multi(:,r), '-o', 'DisplayName', ...
         sprintf('z \\in [%.2f, %.2f]', z_ranges(r,1), z_ranges(r,2)));
end
set(gca, 'YScale', 'log');
xlabel('$\log_{10}(N_{CIV})$', 'Interpreter', 'latex')
ylabel('$\frac{d\mathcal{N}}{d\mathcal{X}_{eff}}$', 'Interpreter','latex')
title('CDDF in Multiple Redshift Ranges');
legend('Location','best');
exportgraphics(fig, 'figSoftP95snr4/phi1D-multi.png', 'Resolution', 400);

