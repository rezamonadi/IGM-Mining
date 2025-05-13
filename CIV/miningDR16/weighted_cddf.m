%% weighted_cddf.m
% Compute weighted C IV Column Density Distribution Function (CDDF)

clear all; close all;

%% 1) Load data
% Parameter catalog
load('processed_qsos_dr16_parameters.mat');

% Likelihood HDF5 via matfile
m      = matfile('processed_qsos_dr16_Likelihood.mat');
h5file = m.Properties.Source;

%% 2) Cosmological setup
Omega_M      = 0.3;
Omega_Lambda = 0.7;
c           = 299792.458;  % km/s
aH0          = 70;         % km/s/Mpc

% Convert to cgs for Omega_CIV later
c_cgs       = c * 1e5;                                 % cm/s
H0_s        = aH0 * 1e5 / 3.0857e24;                   % s^-1
G           = 6.67430e-8;                              % cm^3 g^-1 s^-2
rho_crit    = 3*H0_s^2 / (8*pi*G);                     % g/cm^3

% dX/dz function
dX_dz = @(z) (1 + z).^2 ./ sqrt(Omega_M * (1 + z).^3 + Omega_Lambda);

%% 3) Compute total absorption path length ΔX
min_z = savingCat.all_min_z_c4s(:);
max_z = savingCat.all_max_z_c4s(:);
n_qsos = numel(min_z);

delta_X = zeros(n_qsos, 1);
for i = 1:n_qsos
    if min_z(i) < max_z(i)
        z = linspace(min_z(i), max_z(i), 100);
        delta_X(i) = trapz(z, dX_dz(z));
    end
end

total_dX = sum(delta_X);

%% 4) Set up log N bins
logN_edges   = 12:0.025:17;
dlogN        = logN_edges(2) - logN_edges(1);
logN_centers = logN_edges(1:end-1) + dlogN/2;
nBins        = numel(logN_centers);

%% 5) Read posterior grid & prepare weighted counts
logN_grid = h5read(h5file, '/LikelihoodCat/logN_grid_c4L2');
nGrid     = numel(logN_grid);
info      = h5info(h5file, '/LikelihoodCat/all_sample_log_likelihoods_c4L2');
dims      = info.Dataspace.Size;    % [nGrid, n_qsos, maxAbs]
maxAbs    = dims(3);
n_qsos_lik = dims(2);

weighted_cnt = zeros(1, nBins);

%% 6) Identify valid absorbers (e.g. p_c4 > 0.95)
logN_vec = savingCat.all_map_N_c4L2(:);
p_c4_vec = savingCat.all_p_c4(:);
threshold = 0.95;
valid_lin = find(~isnan(logN_vec) & logN_vec > 0 & p_c4_vec > threshold);
[abs_idx, qso_idx] = ind2sub([maxAbs, n_qsos_lik], valid_lin);

%% 7) Loop over absorbers and accumulate weighted bin counts
for k = 1:numel(valid_lin)
    a = abs_idx(k);
    q = qso_idx(k);

    % load log-likelihood vector for this absorber
    logL = h5read(h5file, '/LikelihoodCat/all_sample_log_likelihoods_c4L2', [1, q, a], [nGrid, 1, 1]);

    % convert to posterior probabilities
    w = exp(logL - max(logL));
    w = w / sum(w);

    % bin the posterior weights
    [~, binID] = histc(logN_grid, logN_edges);
    for ii = 1:nGrid
        b = binID(ii);
        if b >= 1 && b <= nBins
            weighted_cnt(b) = weighted_cnt(b) + w(ii);
        end
    end
end

%% 8) Normalize to get f(N, X)
f_logN_X = weighted_cnt ./ (dlogN * total_dX);
N_centers = 10.^logN_centers;
f_N_X     = f_logN_X ./ (N_centers * log(10));

%% 9) Plot CDDF
figure;
plot(log10(N_centers), log10(f_N_X), 'k+', 'MarkerFaceColor', 'none');
hold on;
xlabel('log_{10}(N_{C IV}/cm^{-2})');
ylabel('log_{10}[f(N,X)]');
title(sprintf('Weighted CDDF (p\_c4 > %.2f)', threshold));
grid on;

%% 10) Fit slope for logN >= 14
fit_min = 14;
valid_fit = isfinite(log10(f_N_X)) & log10(N_centers) >= fit_min;
x_fit = log10(N_centers(valid_fit));
y_fit = log10(f_N_X(valid_fit));
p = polyfit(x_fit, y_fit, 1);
beta = -p(1);
intercept = p(2);

% overlay fit line
y_line = polyval(p, log10(N_centers));
plot(log10(N_centers), y_line, 'r--', 'LineWidth', 1.5);
legend('Data', sprintf('Fit (logN \ge %.2f): slope = %.2f', fit_min, -p(1)));

fprintf('Slope of log-log (logN >= %.2f): %.3f\n', fit_min, -p(1));

%% 11) Compute Omega_CIV from fit
Nmin_term = 10^fit_min;
Nmax_term = max(N_centers);

integral_term = (10^intercept) / (2 - beta) * (Nmax_term^(2-beta) - Nmin_term^(2-beta));
Omega_CIV_fit = (H0_s * (12*1.6726219e-24)) / (c_cgs * rho_crit) * integral_term;

fprintf('Omega_CIV (from fit) = %.3e\n', Omega_CIV_fit);
