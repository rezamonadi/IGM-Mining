% Plot the low-rank GP covariance learned for the DR7 no-C IV QSO model.
%
% Expected input (in the project root by default):
%   learned_model-C13-full.mat
% Output:
%   learned_covariance_dr7.png

model_file = fullfile(pwd, 'learned_model-C13-full.mat');
output_file = fullfile(pwd, 'learned_covariance_dr7.png');

if ~isfile(model_file)
    error('Model file not found: %s', model_file);
end

model = load(model_file, 'M', 'rest_wavelengths');

% M is [number of wavelength pixels x GP rank].  The learned low-rank
% covariance contribution is M*M'.
covariance = model.M * model.M.';

figure('Color', 'w', 'Position', [100 100 900 760]);
imagesc(model.rest_wavelengths, model.rest_wavelengths, covariance);
axis image;
axis xy;
colormap(parula);
colorbar;
xlabel('Rest wavelength (\AA)');
ylabel('Rest wavelength (\AA)');
title(sprintf('DR7 learned GP covariance (rank %d)', size(model.M, 2)));

exportgraphics(gcf, output_file, 'Resolution', 300);
fprintf('Saved covariance plot to: %s\n', output_file);
