%% test_voigt_iP_transition
% Numerical and visual checks for the corrected voigt_iP_transition MEX,
% including a direct comparison with the legacy voigt_iP LSF convention.
%
% Compile first in MATLAB, for example:
mex voigt_iP_transition.c -lcerf -I/path/to/cerf/include -L/path/to/cerf/lib
mex voigt_iP.c -lcerf -I/path/to/cerf/include -L/path/to/cerf/lib

clear;
close all;

if exist('voigt_iP_transition', 'file') == 0
    error(['voigt_iP_transition MEX was not found. Compile it in MATLAB ', ...
           'and ensure its output directory is on the MATLAB path.']);
end
if exist('voigt_iP', 'file') == 0
    error(['voigt_iP MEX was not found. Compile it in MATLAB ', ...
           'and ensure its output directory is on the MATLAB path.']);
end

%% A simple, uniformly sampled fine grid
z = 2.0;
N = 2.0e12;                 % cm^-2
sigma = 5e5;              % cm/s (15 km/s)
fine_factor = 20;            % fine samples per native SDSS pixel
lsf_sigma = 1.10;           % WDISP in native-pixel units
padding = 3 * fine_factor;  % required at both ends by the MEX function
lambda_1548 = 1548.1949462890625;
lambda_1550 = 1550.77001953125;

rest_core = (1543:0.002:1556)';
lambda_core = rest_core * (1 + z);
d_lambda = lambda_core(2) - lambda_core(1);
left_pad = lambda_core(1) - d_lambda * (padding:-1:1)';
right_pad = lambda_core(end) + d_lambda * (1:padding)';
wavelengths = [left_pad; lambda_core; right_pad];
wd_nan = nan(size(wavelengths));       % uses the documented no-LSF path
wd_lsf = lsf_sigma * ones(size(wavelengths));

%% 1. Intrinsic profiles: test all three transition modes
p1548_raw = profile(wavelengths, z, N, 1, sigma, wd_nan, fine_factor);
p1550_raw = profile(wavelengths, z, N, 2, sigma, wd_nan, fine_factor);
pdoublet_raw = profile(wavelengths, z, N, 3, sigma, wd_nan, fine_factor);

assert(numel(p1548_raw) == numel(lambda_core), 'Unexpected output length.');
assert(all(isfinite(pdoublet_raw)), 'Intrinsic doublet contains NaN or Inf.');
assert(all(pdoublet_raw > 0 & pdoublet_raw <= 1), ...
    'An absorption profile must lie in (0, 1].');

% Before convolution, optical depths add, so fluxes multiply exactly.
raw_identity_error = max(abs(pdoublet_raw - p1548_raw .* p1550_raw));
assert(raw_identity_error < 1e-12, ...
    'Mode 3 is not the expected combination of modes 1 and 2.');

[~, i1548] = min(p1548_raw);
[~, i1550] = min(p1550_raw);
assert(abs(rest_core(i1548) - lambda_1548) < 0.03, ...
    'Mode 1 minimum is not near the 1548 transition.');
assert(abs(rest_core(i1550) - lambda_1550) < 0.03, ...
    'Mode 2 minimum is not near the 1550 transition.');

%% 2. Broadened profiles: test the Gaussian LSF branch
p1548_lsf = profile(wavelengths, z, N, 1, sigma, wd_lsf, fine_factor);
p1550_lsf = profile(wavelengths, z, N, 2, sigma, wd_lsf, fine_factor);
pdoublet_lsf = profile(wavelengths, z, N, 3, sigma, wd_lsf, fine_factor);

assert(all(isfinite(pdoublet_lsf)), 'Broadened doublet contains NaN or Inf.');
assert(all(pdoublet_lsf > 0 & pdoublet_lsf <= 1), ...
    'Broadened absorption profile must lie in (0, 1].');
assert(all(pdoublet_lsf <= p1548_lsf + 1e-12), ...
    'Adding the 1550 line should not increase transmitted flux.');
assert(all(pdoublet_lsf <= p1550_lsf + 1e-12), ...
    'Adding the 1548 line should not increase transmitted flux.');

% This is expected to be nonzero: convolution occurs after combining lines.
broadened_difference = max(abs(pdoublet_lsf - p1548_lsf .* p1550_lsf));

%% 3. Regression check: corrected versus legacy oversampled LSF
% The legacy MEX treats WDISP as a fine-grid width and uses only +/- 3
% fine samples.  Its input therefore needs only three fine samples of
% padding, whereas the corrected MEX uses +/- 3 native pixels, or
% +/- 3 * fine_factor fine samples.  Both outputs are aligned to
% lambda_core for this comparison.
legacy_padding = 3;
legacy_wavelengths = [ ...
    lambda_core(1) - d_lambda * (legacy_padding:-1:1)'; ...
    lambda_core; ...
    lambda_core(end) + d_lambda * (1:legacy_padding)' ...
    ];
legacy_wdisp = lsf_sigma * ones(size(legacy_wavelengths));

% Compare the isolated 1548 transition: legacy num_lines = 1 is 1548,
% matching corrected transition_mode = 1.
p1548_legacy = legacy_profile(legacy_wavelengths, z, N, 1, sigma, legacy_wdisp);
assert(numel(p1548_legacy) == numel(lambda_core), ...
    'Legacy and corrected profiles are not aligned to lambda_core.');

line_window = abs(rest_core - lambda_1548) < 0.5;
line_offset = rest_core(line_window) - lambda_1548;
legacy_absorption = 1 - p1548_legacy(line_window);
corrected_absorption = 1 - p1548_lsf(line_window);
legacy_width = sqrt(sum(line_offset.^2 .* legacy_absorption) / sum(legacy_absorption));
corrected_width = sqrt(sum(line_offset.^2 .* corrected_absorption) / sum(corrected_absorption));
legacy_difference = max(abs(p1548_lsf - p1548_legacy));

% With WDISP = 1.10 and fine_factor = 20, the corrected kernel spans
% roughly 22 fine samples per sigma, rather than the legacy 1.10 samples.
% It must consequently be broader and shallower at the isolated line core.
assert(corrected_width > 1.5 * legacy_width, ...
    'Corrected 1548 LSF is not substantially broader than the legacy LSF.');
assert(min(p1548_lsf(line_window)) > min(p1548_legacy(line_window)), ...
    'Corrected 1548 LSF is not shallower than the legacy LSF.');

%% 4. Error paths explicitly implemented by the corrected MEX
must_error(@() profile(wavelengths, z, N, 0, sigma, wd_lsf, fine_factor), ...
    'transition_mode = 0');
must_error(@() profile(wavelengths, z, N, 4, sigma, wd_lsf, fine_factor), ...
    'transition_mode = 4');
must_error(@() profile(wavelengths, z, N, 1.5, sigma, wd_lsf, fine_factor), ...
    'non-integer transition_mode');
must_error(@() profile(wavelengths, z, N, 1, sigma, wd_lsf, 0), ...
    'fine_factor = 0');
must_error(@() profile(wavelengths, z, N, 1, sigma, wd_lsf(1:end-1), fine_factor), ...
    'wrong WDISP length');

fprintf('\nAll numerical assertions passed.\n');
fprintf('Maximum raw doublet identity error: %.3g\n', raw_identity_error);
fprintf('Maximum post-convolution product difference: %.3g\n', broadened_difference);
fprintf('Maximum corrected-versus-legacy 1548 difference: %.3g\n', legacy_difference);
fprintf('Absorption-weighted 1548 width: corrected %.4g A, legacy %.4g A (%.3fx)\n', ...
    corrected_width, legacy_width, corrected_width / legacy_width);

%% Figure 1: each transition mode before instrumental broadening
figure('Name', 'C IV transition modes: intrinsic', 'Color', 'w');
plot(rest_core, p1548_raw, 'LineWidth',1.5, 'Color',[1,0,0,0.7]); hold on;
plot(rest_core, p1550_raw, 'LineWidth', 1.5, 'Color',[0,0,1,0.7]);
plot(rest_core, pdoublet_raw, 'k', 'LineWidth', 1.7, 'Color',[0,0,0,0.7]);
xline(lambda_1548, '--', '1548.1949');
xline(lambda_1550, '--', '1550.7700');
grid on; box on;
xlabel('Rest wavelength (Angstrom)'); ylabel('Transmitted flux');
title('Intrinsic C IV profiles');
legend('mode 1: 1548 only', 'mode 2: 1550 only', 'mode 3: doublet', ...
    'Location', 'best');

%% Figure 2: the effect of instrumental broadening on the doublet
figure('Name', 'C IV doublet: LSF effect', 'Color', 'w');
plot(rest_core, pdoublet_raw, 'k--', 'LineWidth', 1.4); hold on;
plot(rest_core, pdoublet_lsf, 'r', 'LineWidth', 1.7);
grid on; box on;
xlabel('Rest wavelength (Angstrom)'); ylabel('Transmitted flux');
title(sprintf('Instrumental broadening (WDISP = %.2f native pixels)', lsf_sigma));
legend('no LSF (NaN WDISP fallback)', 'Gaussian LSF', 'Location', 'best');

%% Figure 3: zooms confirm that each single-line mode is centered correctly
figure('Name', 'C IV line-center checks', 'Color', 'w');
subplot(1, 2, 1);
plot(rest_core, p1548_raw, 'b', 'LineWidth', 1.5); hold on;
xline(lambda_1548, 'k--'); grid on; box on;
xlim(lambda_1548 + [-0.6, 0.6]); ylim([min(p1548_raw)-0.03, 1.01]);
title('Mode 1: 1548 only'); xlabel('Rest wavelength (Angstrom)'); ylabel('Flux');
subplot(1, 2, 2);
plot(rest_core, p1550_raw, 'm', 'LineWidth', 1.5); hold on;
xline(lambda_1550, 'k--'); grid on; box on;
xlim(lambda_1550 + [-0.6, 0.6]); ylim([min(p1550_raw)-0.03, 1.01]);
title('Mode 2: 1550 only'); xlabel('Rest wavelength (Angstrom)'); ylabel('Flux');

%% Figure 4: demonstrate the order of doublet combination and convolution
figure('Name', 'Convolution order check', 'Color', 'w');
subplot(2, 1, 1);
plot(rest_core, pdoublet_lsf, 'k', 'LineWidth', 1.5); hold on;
plot(rest_core, p1548_lsf .* p1550_lsf, '--', 'LineWidth', 1.2);
grid on; box on; ylabel('Transmitted flux');
legend('broaden(mode 3)', 'broaden(mode 1) .* broaden(mode 2)', ...
    'Location', 'best');
title('Convolution is applied after the doublet is combined');
subplot(2, 1, 2);
plot(rest_core, pdoublet_lsf - p1548_lsf .* p1550_lsf, 'LineWidth', 1.3);
grid on; box on; xlabel('Rest wavelength (Angstrom)'); ylabel('Difference');

%% Figure 5: the corrected LSF must be broader than the legacy fine-grid LSF
figure('Name', 'Corrected versus legacy 1548 LSF', 'Color', 'w');
subplot(2, 1, 1);
plot(rest_core, p1548_legacy, '--', 'LineWidth', 1.4); hold on;
plot(rest_core, p1548_lsf, 'LineWidth', 1.7);
xline(lambda_1548, 'k:'); grid on; box on;
xlim(lambda_1548 + [-0.35, 0.35]);
ylabel('Transmitted flux');
legend('legacy: WDISP used on fine grid', ...
    'corrected: WDISP scaled by fine factor', 'Location', 'best');
title(sprintf('1548 LSF width: corrected / legacy = %.3f', ...
    corrected_width / legacy_width));
subplot(2, 1, 2);
plot(rest_core, p1548_lsf - p1548_legacy, 'LineWidth', 1.3);
xline(lambda_1548, 'k:'); grid on; box on;
xlim(lambda_1548 + [-0.35, 0.35]);
xlabel('Rest wavelength (Angstrom)'); ylabel('Corrected - legacy flux');

%% Local helper functions
function p = profile(wavelengths, z, N, mode, sigma, wdisp, fine_factor)
    p = voigt_iP_transition(wavelengths, z, N, mode, sigma, wdisp, fine_factor);
    p = p(:);
end

function p = legacy_profile(wavelengths, z, N, num_lines, sigma, wdisp)
    p = voigt_iP(wavelengths, z, N, num_lines, sigma, wdisp);
    p = p(:);
end

function must_error(f, description)
    did_error = false;
    try
        f();
    catch
        did_error = true;
    end
    assert(did_error, 'Expected an error for %s.', description);
end
