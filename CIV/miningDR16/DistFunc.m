clear all
close all

% Load the data
load('processed_qsos_dr16_parameters.mat');

% Cosmological parameters
Omega_M = 0.3;
Omega_Lambda = 0.7;

% Constants
c = 299792.458; % km/s
H0 = 70; % Hubble constant in km/s/Mpc

% Define binning for log N (column density)
logN_edges = 12:0.2:17;
logN_centers = logN_edges(1:end-1) + 0.1;
dlogN = 0.2;

% Compute total redshift absorption path ΔX
min_z = savingCat.all_min_z_c4s(:); % Ensure column vector
max_z = savingCat.all_max_z_c4s(:);
n_qsos = length(min_z);

% Function to compute dX/dz
dX_dz = @(z) (1 + z).^2 ./ sqrt(Omega_M * (1 + z).^3 + Omega_Lambda);

% Compute ΔX for each spectrum using trapezoidal rule
delta_X = zeros(n_qsos, 1);
for i = 1:n_qsos
    if min_z(i) < max_z(i)
        z = linspace(min_z(i), max_z(i), 100);
        integrand = dX_dz(z);
        delta_X(i) = trapz(z, integrand);
    else
        delta_X(i) = 0;
    end
end

% Total redshift absorption path
total_dX = sum(delta_X);

% Reshape and extract column density and redshift arrays
logN_all = savingCat.all_map_N_c4L2(:);    % [1D column vector]
z_CIV_all = savingCat.all_map_z_c4L2(:);   % [1D column vector]
p_c4_all = savingCat.all_p_c4(:);
% Construct valid index mask
valid_idx = ~isnan(logN_all) & logN_all > 0 & z_CIV_all > 0 & p_c4_all>0.99;
dirty_idx = ~isnan(logN_all) & logN_all > 0 & z_CIV_all > 0;
% Apply filtering
logN_clean = logN_all(valid_idx);
logN_dirty = logN_all(dirty_idx);
z_CIV_clean = z_CIV_all(valid_idx);

% Bin logN values for f(N, X)
counts = histcounts(logN_clean, logN_edges);
counts_dirty = histcounts(logN_dirty, logN_edges);
% Compute f(N, X)
f_N_X = counts ./ (dlogN * total_dX);
f_n_x_dirty = counts_dirty ./(dlogN* total_dX);
% Plot f(N, X)
figure;
plot(logN_centers, log10(f_N_X), 'ko-', 'MarkerFaceColor', 'k');
xlabel('log_{10}(N_{C IV} / cm^{-2})');
ylabel('f(N, X)');
title('Column Density Distribution Function');
grid on;
hold on 
plot(logN_centers, log10(f_n_x_dirty), 'g.')
% % Plot logN vs z_CIV scatter
% figure;
% scatter(z_CIV_clean, logN_clean, 10, 'b', 'filled');
% xlabel('z_{C IV}');
% ylabel('log_{10}(N_{C IV} / cm^{-2})');
% title('Column Density vs. Redshift of C IV Absorbers');
% grid on;
% % Define binning
z_edges = 1.5:0.2:5;              % Redshift bins
z_centers = z_edges(1:end-1) + 0.1;

logN_edges = 12:0.5:16.5;         % Coarser binning for clearer plot
logN_centers = logN_edges(1:end-1) + 0.25;
dlogN = 0.5;

% Initialize matrix
f_vs_z = zeros(length(logN_centers), length(z_centers));

% % Loop over redshift bins
% for j = 1:length(z_centers)
%     z_min_bin = z_edges(j);
%     z_max_bin = z_edges(j+1);
% 
%     % Select absorbers in redshift bin
%     in_z = z_CIV_clean >= z_min_bin & z_CIV_clean < z_max_bin;
%     logN_in_bin = logN_clean(in_z);
% 
%     % Compute ΔX for the redshift bin
%     z_range = linspace(z_min_bin, z_max_bin, 100);
%     integrand = dX_dz(z_range);
%     delta_X_bin = trapz(z_range, integrand);
% 
%     % Count absorbers in logN bins
%     counts = histcounts(logN_in_bin, logN_edges);
% 
%     % Compute f(N, z) for this z bin
%     f_vs_z(:, j) = counts ./ (dlogN * delta_X_bin);
% end
% 
% % Plot f(N) vs z_CIV for each logN bin
% figure;
% hold on;
% colors = lines(length(logN_centers));

% for i = 1:length(logN_centers)
%     semilogy(z_centers, f_vs_z(i,:), '-o', 'Color', colors(i,:), ...
%              'DisplayName', sprintf('logN = %.1f–%.1f', logN_edges(i), logN_edges(i+1)));
% end
% 
% xlabel('z_{C IV}');
% ylabel('f(N, z)');
% title('Column Density Distribution f(N) vs Redshift');
% legend('Location','northeastoutside');
% grid on;
