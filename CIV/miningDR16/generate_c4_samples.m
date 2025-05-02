% generate_civ_samples: generates CIV parameter samples from training
% catalog

% since we don't have hash tables in catalog.mat, we load the ascii file directly




c4_catalog = load(sprintf('%s/c4_catalog', c4_catalog_directory('Cooksey_C4_cat')));
c4_catalog =c4_catalog(c4_catalog(:,3)>0,:); % removing some null column densities
% training_set_name = 'UVES';

% c4_catalog = load(sprintf('%s/c4_catalog', c4_catalog_directory(training_set_name)));
% c4_catalog = fitsread('data/C4_catalogs/UVES_C4_cat/tab2.fits', 'binarytable');

% generate quasirandom samples from p(normalized offset, log₁₀(N_CIV))
rng('default');
sequence = scramble(haltonset(3), 'rr2');

% the first dimension can be used directly for the uniform prior over
% offsets
offset_z_samples  = sequence(1:num_C4_samples, 1)';
offset_sigma_samples = sequence(1:num_C4_samples, 1)';
% we must transform the second dimension to have the correct marginal
% distribution for our chosen prior over column density, which is a
% mixture of a uniform distribution on log₁₀ N_CIV and a distribution
% we fit to observed data

% uniform component of column density prior
u = makedist('uniform', ...
             'lower', uniform_min_log_nciv, ...
             'upper', uniform_max_log_nciv);

% extract observed log₁₀ N_CIV samples directly from CIV catalog
log_nciv = c4_catalog(:, 3);
% log_nciv = c4_catalog{6};

% make a quadratic fit to the estimated log p(log₁₀ N_CIV) over the
% specified range
x = linspace(fit_min_log_nciv, fit_max_log_nciv, 1e3);
kde_pdf = ksdensity(log_nciv, x);
f = polyfit(x, log(kde_pdf), 2);

% solve for the turning point of the quadratic function:
% f = a x^2 + b x + c; f' = a * 2 x + b
% below the turning point, we assume a uniform prior.
% This is due to lower column density lines are harder to
% measure, so they are very possible incomplete.
turning_log_nciv = roots([ 2*f(1) f(2) ]);
fprintf('Solved roots( df/dN ) : %d\n', turning_log_nciv);

% convert this to a PDF and normalize
% extrapolate the value at the turning point (~14.4) to 12.5 < logNCIV < 14.4
unnormalized_pdf = ...
     @(nciv) ( exp(polyval(f,  nciv))              .*      heaviside( nciv - turning_log_nciv ) ...
           +   exp(polyval(f,  turning_log_nciv))  .* (1 - heaviside( nciv - turning_log_nciv )) );
Z = integral(unnormalized_pdf, fit_min_log_nciv, 18); 
% integrate until 16 to get the tail region

% create the PDF of the mixture between the unifrom distribution and
% the distribution fit to the data
normalized_pdf = @(nciv) ...
          alpha  * (unnormalized_pdf(nciv) / Z) + ...
     (1 - alpha) * (pdf(u, nciv));

cdf = @(nciv) (integral(normalized_pdf, fit_min_log_nciv, nciv));

% use inverse transform sampling to convert the quasirandom samples on
% [0, 1] to appropriate values
log_nciv_samples = zeros(1, num_C4_samples);
W_r_1548_samples = nan([num_C4_samples,1]);
W_r_1550_samples = nan([num_C4_samples,1]);

f = waitbar(0, 'Starting');
for i = 1:num_C4_samples
     log_nciv_samples(i) = ...
          fzero(@(nciv) (cdf(nciv) - sequence(i, 3)), 14.4);
     waitbar(i/num_C4_samples, f, sprintf('Progress: %d %%',...
      floor(i/num_C4_samples*100)));
end

% precompute N_CIV samples for convenience
nciv_samples = 10.^log_nciv_samples;


variables_to_save = {'offset_z_samples', 'offset_sigma_samples', 'log_nciv_samples', 'nciv_samples'};
             
 save(sprintf('%s/civ_samples_N-%d-%d-Sigma-%d-%d-Num-%d.mat', processed_directory(releaseTest),...
     floor(uniform_min_log_nciv*100),floor(uniform_max_log_nciv*100), floor(min_sigma/1e5), floor(max_sigma/1e5), num_C4_samples),  variables_to_save{:}, '-v7.3');

fig = figure();

histogram(log_nciv_samples)
exportgraphics(fig, sprintf('%s-hist2D.png', training_set_name), 'Resolution', 800)

