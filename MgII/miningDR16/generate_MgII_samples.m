% generate_civ_samples: generates CIV parameter samples from training
% catalog

% since we don't have hash tables in catalog.mat, we load the ascii file directly



MgII_catalog = load('data/dr7/processed/MgII-cat.mat');

% generate quasirandom samples from p(normalized offset, log₁₀(N_CIV))
rng('default');
sequence = scramble(haltonset(3), 'rr2');

% the first dimension can be used directly for the uniform prior over
% offsets
offset_z_samples  = sequence(1:num_MgII_samples, 1)';
offset_sigma_samples = sequence(1:num_MgII_samples, 2)';
% we must transform the second dimension to have the correct marginal
% distribution for our chosen prior over column density, which is a
% mixture of a uniform distribution on log₁₀ N_CIV and a distribution
% we fit to observed data

% uniform component of column density prior
u = makedist('uniform', ...
             'lower', uniform_min_log_nMgII, ...
             'upper', uniform_max_log_nMgII);

% extract observed log₁₀ N_CIV samples directly from CIV catalog
log_nMgII = MgII_catalog.NMgII;
% log_nMgII = c4_catalog{6};

% make a quadratic fit to the estimated log p(log₁₀ N_CIV) over the
% specified range
x = linspace(fit_min_log_nMgII, fit_max_log_nMgII, 1e3);
kde_pdf = ksdensity(log_nMgII, x);
f = polyfit(x, log(kde_pdf), 2);

% solve for the turning point of the quadratic function:
% f = a x^2 + b x + c; f' = a * 2 x + b
% below the turning point, we assume a uniform prior.
% This is due to lower column density lines are harder to
% measure, so they are very possible incomplete.
turning_log_nMgII = roots([ 2*f(1) f(2) ]);
fprintf('Solved roots( df/dN ) : %d\n', turning_log_nMgII);

% convert this to a PDF and normalize
% extrapolate the value at the turning point (~14.4) to 12.5 < lognMgII < 14.4
unnormalized_pdf = ...
     @(nMgII) ( exp(polyval(f,  nMgII))              .*      heaviside( nMgII - turning_log_nMgII ) ...
           +   exp(polyval(f,  turning_log_nMgII))  .* (1 - heaviside( nMgII - turning_log_nMgII )) );
Z = integral(unnormalized_pdf, fit_min_log_nMgII, 18); 
% integrate until 16 to get the tail region

% create the PDF of the mixture between the unifrom distribution and
% the distribution fit to the data
normalized_pdf = @(nMgII) ...
          alpha  * (unnormalized_pdf(nMgII) / Z) + ...
     (1 - alpha) * (pdf(u, nMgII));

cdf = @(nMgII) (integral(normalized_pdf, fit_min_log_nMgII, nMgII));

% use inverse transform sampling to convert the quasirandom samples on
% [0, 1] to appropriate values
log_nMgII_samples = zeros(1, num_MgII_samples);


f = waitbar(0, 'Starting');
for i = 1:num_MgII_samples
     log_nMgII_samples(i) = ...
          fzero(@(nMgII) (cdf(nMgII) - sequence(i, 3)), 14.4);
     waitbar(i/num_MgII_samples, f, sprintf('Progress: %d %%',...
      floor(i/num_MgII_samples*100)));
end

% precompute N_CIV samples for convenience
nMgII_samples = 10.^log_nMgII_samples;

variables_to_save = {'offset_z_samples', 'offset_sigma_samples', 'log_nMgII_samples', 'nMgII_samples'};
             
 save(sprintf('%s/MgII_samples_%s.mat', processed_directory(releaseTest), sample_name),  variables_to_save{:}, '-v7.3');

fig = figure();

histogram(log_nMgII_samples)
exportgraphics(fig, sprintf('%s-hist2D.png', sample_name), 'Resolution', 800)

