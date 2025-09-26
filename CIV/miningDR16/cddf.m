%---------- C IV CDDF ---------------
% clearvars; close all; 
% clc;
SNR_threshold = 4;
p_threshold = 0.85;
nSamples=1000;
completeness = false; % set to false if you want to load precomputed completeness
loading = false; % set to false if you want to load pre-saved sampledAbsorbers and processed_qsos_dr16_parameters.mat
cap =0.3;
% % loading data

if loading
    load(sprintf('sampled_logN_per_absorber_Smp%d.mat',nSamples));   % contains sampledAbsorbers struct array
    load('processed_qsos_dr16_parameters.mat');      % contains savingCat struct array
    catDR16 = load('catalog.mat');
    if completeness 
        load('preloaded_qsos_DR16.mat');  % contains spec struct array
    end
end

%% 

figDir = sprintf('Figs-snr%d-p%f-smp%d-capG%f', SNR_threshold, p_threshold, nSamples, cap);
if ~isfolder(figDir)
    mkdir(figDir)
end
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

log10_N_edges = linspace(min(logN_samples,[],'all')-0.25, max(logN_samples,[],'all')+0.25, 10);
z_edges       = linspace(min(z_abs_samples,[],'all')-0.25, max(z_abs_samples,[],'all')+0.25, 30);
opts.addPoisson = false;  % add Poisson counting variance on top of soft-count variance
Out = soft_hist_2d(p_c4, logN_samples, z_abs_samples, log10_N_edges, z_edges, opts); 



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


%% Completeness calculation with b-marginalization
% function out = civ_sensitivity_gNz_b1D(all_wavelengths, ...
%                                        all_sigma_pixel, ...
%                                        all_pixel_mask, ...
%                                        all_noise_vars, ...
%                                        N_grid, z_edges, b_grid_kms, b_weights, opts)
if completeness
    wave = all_wavelengths(filter_flags==0);
    px_sigma = all_sigma_pixel(filter_flags==0);
    px_mask  = all_pixel_mask(filter_flags==0);
    noise_var = all_noise_variance(filter_flags==0);
    h = histogram(savingCat.all_map_sigma_c4L2(savingCat.all_p_c4>0.85),'NumBins',50, 'Normalization','probability');
    bWeight = h.Values;
    bGrid = h.BinEdges(1:end-1) + 0.5*h.BinWidth;
    opts = struct;
    OutPut = civ_sensitivity_gNz_b1D(wave, ...
                                    px_sigma, ...
                                    px_mask, ...
                                    noise_var,  ...
                                    10.^(Out.xcenters), ...
                                    z_edges, ...
                                    bGrid, ...
                                    bWeight, ...
                                    opts);

    save('completeness_output.mat', 'OutPut');
else
    load('completeness_output.mat'); % load precomputed completeness
end

%% --- Build the raw 2D counts (log10 N vs z) ---
% (Rebuild to ensure counts match h.BinWidth and edges)
counts     = Out.H2D;  % size = nLogN x nZ
var_counts = (Out.sigma2D).^2; % variance from soft histogram
logN_ctr   = Out.xcenters; % assuming uniform bin width
z_ctr      = Out.ycenters; % assuming uniform bin width
Delta_logN = Out.dx(1);  % assuming uniform bin width
Delta_z    = Out.dy(1);  % assuming uniform bin width

% Convert each (logN,z) bin: φ -> f before averaging
g_counts = OutPut.gNz;  % completeness grid [nN x nZ]

% ----- 1) Convert g_counts -> probability in [0,1] (column-wise) -----
g_prob = zeros(size(g_counts));
nZ = size(g_counts,2);
for k = 1:nZ
    Tk = max(g_counts(:,k));             % plateau / total opportunities in this z-bin
    if Tk > 0
        g_prob(:,k) = g_counts(:,k) / Tk;
    else
        g_prob(:,k) = 0;
    end
end
g_prob = min(max(g_prob,0),1);           % clamp numeric noise to [0,1]
% ----- 2) Low-g handling (cap or mask) -----
g_eff = max(g_prob, cap);         % cap weights via g_eff

C_of_N = OutPut.C_of_N; % completeness per N [nN x 1]
C_of_N = min(max(C_of_N, 0), 1); % clamp numeric noise to [0,1]
C_of_N(C_of_N < cap) = cap;
% C_of_N(C_of_N > 0.9) = 0.9;
% --- Normalize to φ(logN, z) = counts / (ΔlogN * ΔX_eff) ---

phi_logN_z = zeros(size(counts));
phi_logN_z_comp_cor = zeros(size(counts));
var_phi_logN_z = zeros(size(counts));
var_phi_logN_z_com_cor = zeros(size(counts));
for k = 1:numel(z_ctr)
    dX_eff = DeltaX_eff_per_zbin(k);
    if dX_eff > 0 && isfinite(dX_eff)
        denom               = (Delta_logN * dX_eff);
        
        phi_logN_z_comp_cor(:,k)     = counts(:,k)    / denom./ g_eff(:,k);  % <-- divide by completeness
        var_phi_logN_z_com_cor(:,k) = var_counts(:,k) / (denom.^2)./g_eff(:,k).^2;   % <-- square the denom

        phi_logN_z(:,k)     = counts(:,k)    / denom;  % <-- divide by completeness
        var_phi_logN_z(:,k) = var_counts(:,k) / (denom.^2);   % <-- square the denom

    else
        phi_logN_z_comp_cor(:,k) = 0;
        var_phi_logN_z_com_cor(:,k) = 0;

        phi_logN_z(:,k) = 0;
        var_phi_logN_z(:,k) = 0;
    end
end



% % %% --- Plot φ(logN, z) as heatmap ---
% % fig = figure;
% % imagesc(logN_ctr, z_ctr, phi_logN_z.' ./ (log(10)*(10.^logN_ctr.'))); % convert to per dex
% % colorbar;
% % xlabel('$log_{10}(N_{CIV})$', 'Interpreter', 'latex');
% % ylabel('$z_{\rm{CIV}}$', 'Interpreter', 'latex');
% % title('$\phi(log_{10}(N_{CIV}), z)$  [per dlogN per dX]', 'Interpreter','latex');
% % exportgraphics(fig, sprintf('%s/phi.png', fidDir), 'Resolution', 400)

%% --- Redshift-averaged 1D CDDF over your full z-range (path-length weighted) ---
DeltaX = DeltaX_eff_per_zbin(:);   % [nZ x 1]
X_tot  = sum(DeltaX);

Nlin        = 10.^logN_ctr(:);                 % [nN x 1]
% phi_to_f    = 1 ./ (Nlin * log(10));           % [nN x 1]
phi_to_f    = 1 ./ (Nlin*log(10));           % [nN x 1]
phi_to_f_sq = phi_to_f.^2;





f_bins     = phi_logN_z .* phi_to_f;                          % [nN x nZ]
var_f_bins = var_phi_logN_z .* phi_to_f_sq;                   % [nN x nZ]

f_bins_comp_cor    = phi_logN_z_comp_cor .* phi_to_f;                          % [nN x nZ]
var_f_bins_comp_cor = var_phi_logN_z_com_cor .* phi_to_f_sq;                   % [nN x nZ]


% Average over z with ΔX weights
f_1D     = (f_bins * DeltaX) / X_tot;                         % [nN x 1]
err_f_1D_1 = sqrt( (var_f_bins * (DeltaX.^2)) ) / X_tot;        % [nN x 1]
err_f_1D_2 = Out.sigmaHx ./ (Delta_logN * X_tot) ./ (Nlin* log(10));  % Poisson error from total counts, converted to f

f_1D_comp_cor     = (f_bins_comp_cor * DeltaX) / X_tot;                         % [nN x 1]
err_f_1D_1_comp_cor = sqrt( (var_f_bins_comp_cor * (DeltaX.^2)) ) / X_tot;        % [nN x 1]
err_f_1D_2_comp_cor = Out.sigmaHx ./ (Delta_logN * X_tot) ./ (Nlin* log(10))./C_of_N;  % Poisson error from total counts, converted to f


% Plot (per-N)
figure;
errorbar(logN_ctr, f_1D, err_f_1D_1, '-', 'LineWidth', 2); hold on
errorbar(logN_ctr, f_1D_comp_cor, err_f_1D_1_comp_cor, '-', 'LineWidth',2); hold on

xlabel('$\log_{10} N_{\rm C\,IV}$','Interpreter','latex');
ylabel('$f(N)\ [{\rm cm}^{2}\ X^{-1}]$','Interpreter','latex');
grid on
title(sprintf('CDDF (per N): SNR>%d, p>%.2f, %d samples', SNR_threshold, p_threshold, nSamples));
set(gca,'YScale','log');
set(gca, 'FontSize', 20);




% Beta 

logN_fit_min = 13.5; logN_fit_max = 15.5;

I_fit = (logN_ctr >= logN_fit_min) & (logN_ctr <= logN_fit_max) & (f_1D_comp_cor.' > 0) & (err_f_1D_1_comp_cor.' > 0); 

sigma_raw = err_f_1D_2_comp_cor(I_fit);           % the one you pass to the fitter
f_sel     = f_1D_comp_cor(I_fit);

rel_floor = 0.30;                                  % 20% relative floor (tune 0.1–0.3)
nz        = sigma_raw(sigma_raw>0);
abs_floor = 0.5*median(nz);

sigma_eff = sqrt( sigma_raw.^2 + (rel_floor*abs(f_sel)).^2 );
sigma_eff = max(sigma_eff, abs_floor);

[beta,f0,beta_err,f0_err,stats] = fit_cddf_powerlaw(Nlin(I_fit), f_sel, sigma_eff);
N0 = stats.N0;
plot(logN_ctr(I_fit), f0*(Nlin(I_fit)/N0).^(-beta), 'r--', 'LineWidth',2);

legend('Raw', 'Completeness-corrected', 'Fit', 'Location','best');

exportgraphics(gcf, sprintf('%s/CDDF_perN-CompCorrected.png', figDir), 'Resolution', 400);
% % % Omega all with CDDF 

% % % [Omega, sigmaOmega, pieces] = omegaCIV_from_powerlaw(beta, beta_err, f0, f0_err, 1e14, 1e13, 1e14, opts)


% %% line density 
% fig=figure;

% dNdX_per_z = sum(phi_logN_z, 1)*Delta_logN;
% err_dNdX_per_z = Out.sigmaHy./DeltaX_eff_per_zbin.';
% errorbar(z_ctr, (dNdX_per_z), err_dNdX_per_z); yscale('log')
% exportgraphics(fig, sprintf('%s/LineDensity.png', figDir), 'Resolution', 1000)

%% CDDF for custom z-ranges (exact ΔX weighting, no linear-in-z assumption)

% Prereqs:
%   z_edges                % [1 x (nZ+1)]  same edges used to build phi_logN_z
%   nSL, z_min_clipped, z_max_clipped, telluric_z, int_dX
%   phi_logN_z             % [nN x nZ]     per-dex CDDF in each (logN,z) bin
%   var_phi_logN_z         % [nN x nZ]     variance of phi in each bin (else zeros)
%   logN_ctr               % [nN x 1]
%   Delta_logN             % scalar (dex)   (not used here, but handy later)

if ~exist('var_phi_logN_z','var') || isempty(var_phi_logN_z)
    var_phi_logN_z = zeros(size(phi_logN_z));
end

nZ      = numel(z_edges) - 1;
z_ranges = [0.5 1.5; 1.5 2.0; 2.0 3.0; 3.0 4.0; 4.0 5.0];  % EDIT as needed
nRanges = size(z_ranges,1);

phi_1D_multi      = zeros(numel(logN_ctr), nRanges);   % per-dex CDDF
err_phi_1D_multi  = zeros(numel(logN_ctr), nRanges);

% W_range(r,j): ΔX contributed by z-bin j inside z-range r
W_range = zeros(nRanges, nZ);

for r = 1:nRanges
    z_lo_r = z_ranges(r,1);
    z_hi_r = z_ranges(r,2);

    % exact ΔX inside the requested range, per z-bin j
    for j = 1:nZ
        zb_lo = z_edges(j);
        zb_hi = z_edges(j+1);
        lo = max(zb_lo, z_lo_r);
        hi = min(zb_hi, z_hi_r);
        if hi <= lo, continue; end

        dX_sum = 0;
        for s = 1:nSL
            zlo = max(z_min_clipped(s), lo);
            zhi = min(z_max_clipped(s), hi);
            if zhi <= zlo, continue; end
            dX_s = int_dX(zlo, zhi);
            % subtract tellurics per sightline
            for t = 1:size(telluric_z,1)
                zt_lo = max(zlo, telluric_z(t,1));
                zt_hi = min(zhi, telluric_z(t,2));
                if zt_hi > zt_lo
                    dX_s = dX_s - int_dX(zt_lo, zt_hi);
                end
            end
            dX_sum = dX_sum + dX_s;
        end
        W_range(r,j) = dX_sum;
    end

    X_range = sum(W_range(r,:));
    if X_range <= 0
        warning('No path in z-range [%.2f, %.2f].', z_lo_r, z_hi_r);
        continue
    end

    % ---- Per-dex CDDF in this range (path-length–weighted over z) ----
    % phi_range = (phi_logN_z * W_range(r,:).') / X_range;   % [nN x 1]
    weights = (W_range(r,:).');                               % [nZ x 1]
    phi_r   = (phi_logN_z * weights) / X_range;               % [nN x 1]
    phi_1D_multi(:,r) = phi_r;

    % ---- Uncertainty propagation (assuming bin variances add with weights^2) ----
    % Var[Σ w_j φ_j] = Σ w_j^2 Var[φ_j]; then divide by X_range^2.
    var_phi_r = (var_phi_logN_z * (weights.^2)) / (X_range^2);  % [nN x 1]
    err_phi_1D_multi(:,r) = sqrt(var_phi_r);
end

% Plot per-dex CDDFs
fig = figure; hold on; box on; grid on;
for r = 1:nRanges
    y = phi_1D_multi(:,r);
    e = err_phi_1D_multi(:,r);
    if all(y==0), continue; end
    errorbar(logN_ctr, y, e, '-o', 'DisplayName', ...
        sprintf('z \\in [%.2f, %.2f]', z_ranges(r,1), z_ranges(r,2)));
end
set(gca,'YScale','log');
xlabel('$\log_{10} N_{\rm C\,IV}$','Interpreter','latex');
ylabel('$\phi(\log N)\ [{\rm dex^{-1}}\ X^{-1}]$','Interpreter','latex');
title('CDDF in Multiple Redshift Ranges (exact $\Delta X$ weighting)','Interpreter','latex');
legend('Location','best');
exportgraphics(fig, sprintf('%s/CDDF-Multi.png', figDir), 'Resolution', 400);


[Omega_z, sigmaOmega_z, z_mid, out] = omegaCIV_slices_exactDX( ...
    phi_logN_z, var_phi_logN_z, logN_ctr, Delta_logN, ...
    z_edges, z_ranges, ...
    nSL, z_min_clipped, z_max_clipped, telluric_z, int_dX);

fig = figure; hold on; box on; grid on;
errorbar(z_mid, Omega_z, sigmaOmega_z, '-o');
xlabel('$z$','Interpreter','latex');
ylabel('$\Omega_{\rm C\,IV}\ [10^{-8}]$','Interpreter','latex');
title('Cosmic Mass Density of C IV','Interpreter','latex');
exportgraphics(fig, sprintf('%s/OmegaCIV-Multi.png', figDir), 'Resolution', 400);