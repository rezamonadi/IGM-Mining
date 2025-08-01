% learn_qso_model: fits GP to training catalog via maximum likelihood

rng('default');

% load catalog
catalog = load(sprintf('%s/catalog', processed_directory(release)));

% load preprocessed QSOs
variables_to_load = {'all_wavelengths', 'all_flux', 'all_noise_variance', ...
                     'all_pixel_mask', 'all_num_c4'};
preqsos = matfile(preloaded_qsos_cat_name);

% determine which spectra to use for training; allow string value for train_ind
if learning ~= 2
  if (ischar(train_ind))
    train_ind = eval(train_ind);
  end
end

% select training vectors
all_wavelengths    = preqsos.all_wavelengths;
if learning == 2 
  train_ind = true(size(all_wavelengths));
end

all_wavelengths    = all_wavelengths(train_ind, :);
all_flux           = preqsos.all_flux;
all_flux           = all_flux(train_ind, :);
all_noise_variance = preqsos.all_noise_variance;
all_noise_variance = all_noise_variance(train_ind, :);
all_pixel_mask     = preqsos.all_pixel_mask;
all_pixel_mask     = all_pixel_mask(train_ind, :);
z_qsos             = catalog.all_zqso(train_ind);

clear preqsos

num_quasars = numel(z_qsos);
rest_wavelengths = (min_lambda:dlambda:max_lambda);
num_rest_pixels  = numel(rest_wavelengths);

rest_fluxes          = nan(num_quasars, num_rest_pixels);
rest_noise_variances = nan(num_quasars, num_rest_pixels);
is_empty             = false(num_quasars, 1);

norm_window = (rest_wavelengths > 2200 & rest_wavelengths < 2300);

for i = 1:num_quasars
  z_qso = z_qsos(i);

  this_wavelengths    = all_wavelengths{i}';
  this_flux           = all_flux{i}';
  this_noise_variance = all_noise_variance{i}';
  this_pixel_mask     = all_pixel_mask{i}';

  this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qso);

  this_flux(this_pixel_mask) = nan;
  this_noise_variance(this_pixel_mask) = nan;

  if isempty(this_wavelengths)
    is_empty(i) = true;
    continue;
  end

  interp_flux = interp1(this_rest_wavelengths, this_flux, rest_wavelengths);
  interp_noise = interp1(this_rest_wavelengths, this_noise_variance, rest_wavelengths);

  median_flux = nanmedian(interp_flux(norm_window));
  if isnan(median_flux) || median_flux <= 0
    is_empty(i) = true;
    continue;
  end

  interp_flux = interp_flux / median_flux;
  interp_noise = interp_noise / (median_flux^2);

  rest_fluxes(i, :) = interp_flux;
  rest_noise_variances(i, :) = interp_noise;

  fprintf('Processed quasar %i/%i\n', i, num_quasars);
end

% filter out empty spectra
z_qsos               = z_qsos(~is_empty);
rest_fluxes          = rest_fluxes(~is_empty, :);
rest_noise_variances = rest_noise_variances(~is_empty, :);
num_quasars = numel(z_qsos);

fprintf('Filtered out %i empty spectra. Remaining quasars: %i\n', sum(is_empty), num_quasars);

% mask noisy pixels
ind = (rest_noise_variances > max_noise_variance);
fprintf("Masking %g of pixels\n", nnz(ind) / numel(ind));
rest_fluxes(ind) = nan;
rest_noise_variances(ind) = nan;

% Filter out spectra with too many NaNs
ind = sum(isnan(rest_fluxes), 2) < num_rest_pixels - min_num_pixels;
fprintf("Filtering %g quasars for NaNs\n", num_quasars - nnz(ind));
rest_fluxes = rest_fluxes(ind, :);
rest_noise_variances = rest_noise_variances(ind, :);

% Check for wavelength bins with too many NaNs
nancolfrac = sum(isnan(rest_fluxes), 1) / nnz(ind);
fprintf("Columns with nan > 0.9: ");
max(find(nancolfrac > 0.9))

% find empirical mean and center data
mu = nanmean(rest_fluxes);
centered_rest_fluxes = bsxfun(@minus, rest_fluxes, mu);
clear('rest_fluxes');

% Fill NaNs for PCA only
pca_centered_rest_flux = centered_rest_fluxes;
[num_quasars, ~] = size(pca_centered_rest_flux);
for i = 1:num_quasars
  row = pca_centered_rest_flux(i, :);
  row(isnan(row)) = nanmedian(row);
  pca_centered_rest_flux(i, :) = row;
end

% PCA initialization
[coefficients, ~, latent] = pca(pca_centered_rest_flux, 'numcomponents', k, 'rows', 'complete');
initial_M = bsxfun(@times, coefficients(:, 1:k), sqrt(latent(1:k))');

% GP model fitting via L-BFGS
objective_function = @(x) objective(x, centered_rest_fluxes, rest_noise_variances);
[x, log_likelihood, ~, minFunc_output] = minFunc(objective_function, initial_M, minFunc_options);

% Reshape learned basis
ind = (1:(num_rest_pixels * k));
M = reshape(x(ind), [num_rest_pixels, k]);

% Save model
if learning == 1
  variables_to_save = {'release', 'train_ind', 'max_noise_variance', ...
                      'minFunc_options', 'rest_wavelengths', 'mu', ...
                      'initial_M', 'M',  'log_likelihood', ...
                      'minFunc_output', 'test_ind', 'prior_ind'};
  save(sprintf('%s/learned_model-%s.mat', processed_directory(release), training_set_name), ...
       variables_to_save{:}, '-v7.3');
end

if learning == 2
  variables_to_save = {'release', 'max_noise_variance', 'minFunc_options', ...
                       'rest_wavelengths', 'mu', 'initial_M', 'M', ...
                       'log_likelihood', 'minFunc_output'};
  save(sprintf('%s/learned_model-C13-full.mat', processed_directory(release)), ...
       variables_to_save{:}, '-v7.3');
end
