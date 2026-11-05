function report = diagnose_civ_single_run(result_file, varargin)
%DIAGNOSE_CIV_SINGLE_RUN Validate one saved DR16 C IV processing run.
%
% This is deliberately a single-file diagnostic: it does not compare to a
% second nAVG run and does not require the original proposal MAT-file.  It is
% intended for outputs such as the 5,000-QSO run, whose saved compact posterior
% draws replace the full REW sample cube used by diagnose_civ_pilot.
%
% Example:
%   report = diagnose_civ_single_run( ...
%       'processed_qsos_tst_DR16-8_S_1_E_5000.mat', 'MakePlots', true);
%
% The full marginal-likelihood check reads every saved L2 likelihood sample.
% Set 'CheckSampleLikelihoods', false only for a quick structural check.
% With MakePlots=true, a second figure pools the compact posterior draws from
% cells passing each p_c4 threshold.  Each posterior quantity has one panel;
% threshold curves within that panel use identical histogram bins.

if nargin < 1 || isempty(result_file)
    result_file = 'processed_qsos_tst_DR16-8_S_1_E_5000.mat';
end
p = inputParser;
p.addParameter('BlockQSOs', 100, @(x) isscalar(x) && x >= 1 && x == floor(x));
p.addParameter('PC4Threshold', [0.50 0.75 0.85 0.95 0.99], ...
    @(x) isnumeric(x) && isvector(x) && all(isfinite(x)) && all(x >= 0 & x <= 1));
p.addParameter('CheckSampleLikelihoods', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MakePlots', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

assert(isfile(result_file), 'Result file not found: %s', result_file);
m = matfile(result_file);
names = who(m);
required = {'p_no_c4', 'p_c4L1', 'p_c4', ...
    'log_priors_no_c4', 'log_priors_c4', ...
    'log_likelihoods_no_c4', 'log_likelihoods_c4L1', 'log_likelihoods_c4L2', ...
    'log_posteriors_no_c4', 'log_posteriors_c4L1', 'log_posteriors_c4L2', ...
    'map_z_c4L2', 'map_N_c4L2', 'map_sigma_c4L2', ...
    'min_z_c4s', 'max_z_c4s', ...
    'posterior_z_c4L2_draws', 'posterior_logN_c4L2_draws', ...
    'posterior_sigma_c4L2_draws', 'posterior_REW_1548_draws', ...
    'posterior_REW_1550_draws', 'posterior_z_c4L2_16', 'posterior_z_c4L2_84', ...
    'posterior_logN_c4L2_16', 'posterior_logN_c4L2_84', ...
    'posterior_sigma_c4L2_16', 'posterior_sigma_c4L2_84', ...
    'posterior_REW_1548_16', 'posterior_REW_1548_84', ...
    'posterior_REW_1550_16', 'posterior_REW_1550_84'};
missing = required(~ismember(required, names));
assert(isempty(missing), 'Required saved array(s) missing: %s', strjoin(missing, ', '));

shape = size(m, 'p_c4');
assert(numel(shape) == 2, 'p_c4 must be a [QSO, component] array.');
nq = shape(1); nc = shape(2);
assert(nq > 0 && nc > 0, 'p_c4 is empty.');
check_shapes(m, required, nq, nc);
ndraw = size(m, 'posterior_z_c4L2_draws'); ndraw = ndraw(2);

fprintf('\nSingle-run C IV diagnostic: %s\n', result_file);
fprintf('Saved output: %d QSOs x %d components; %d posterior draws/cell.\n', nq, nc, ndraw);

% These small arrays are read once.  A cell is active only when all three
% model probabilities are finite; partial cells are always a corruption sign.
p0 = m.p_no_c4; p1 = m.p_c4L1; p2 = m.p_c4;
active_any = isfinite(p0) | isfinite(p1) | isfinite(p2);
active = isfinite(p0) & isfinite(p1) & isfinite(p2);
probability_invalid = nnz(active_any & ~active) + nnz(active & ...
    (p0 < -1e-12 | p1 < -1e-12 | p2 < -1e-12 | p0 > 1+1e-12 | ...
     p1 > 1+1e-12 | p2 > 1+1e-12 | abs(p0 + p1 + p2 - 1) > 1e-10));

% Searches are iterative: active component cells must form a prefix per QSO.
prefix_violation = false(nq, 1);
for q = 1:nq
    a = find(active(q,:), 1, 'last');
    if ~isempty(a), prefix_violation(q) = any(~active(q,1:a)); end
end
active_by_component = sum(active, 1);
processed_qsos = sum(any(active, 2));
fprintf('Active probability cells: %d; QSOs with >=1 completed component: %d/%d.\n', ...
    nnz(active), processed_qsos, nq);
fprintf('Completed QSOs by component: %s\n', mat2str(active_by_component));

algebra = check_posterior_algebra(m, active, nq, nc);
map = check_map_and_intervals(m, active, nq, nc, ndraw, opt.BlockQSOs);
rew = check_rew_draws(m, p2, active, nq, nc, opt.PC4Threshold, opt.BlockQSOs, opt.MakePlots);
histograms = struct('made', false);
if opt.MakePlots
    histograms = plot_threshold_posteriors(m, p2, active, nq, nc, opt.PC4Threshold, opt.BlockQSOs);
end
likelihood = check_sample_likelihoods(m, names, active, nq, nc, opt);
thresholds = threshold_table(p2, active, opt.PC4Threshold);

report = struct;
report.result_file = result_file;
report.shape = [nq nc];
report.posterior_draws = ndraw;
report.missing_arrays = missing;
report.active_cells = nnz(active);
report.processed_qsos = processed_qsos;
report.active_by_component = active_by_component;
report.probability_invalid = probability_invalid;
report.prefix_violation_qsos = find(prefix_violation);
report.posterior_algebra = algebra;
report.map_and_intervals = map;
report.rew = rew;
report.posterior_histograms = histograms;
report.sample_likelihoods = likelihood;
report.p_c4_thresholds = thresholds;

numerical_failures = probability_invalid + nnz(prefix_violation) + algebra.invalid + ...
    map.invalid + likelihood.mismatch;
report.pass_numerical = (numerical_failures == 0);
% A physical REW-ratio flag is reported separately: it can reveal a modeling
% issue, but it is not evidence that a v7.3 write was interrupted or corrupt.
report.pass_rew_ratio = (rew.negative == 0 && rew.stronger1550 == 0 && rew.ratio_outside == 0);

fprintf('\nProbability checks: %d invalid/partial cells.\n', probability_invalid);
fprintf('Posterior algebra: %s\n', algebra.message);
fprintf('MAP and intervals: %s\n', map.message);
fprintf('Marginal likelihoods: %s\n', likelihood.message);
fprintf('Posterior REWs: finite=%d, negative=%d, W1550>W1548=%d, ratio outside [1,2]=%d.\n', ...
    rew.finite, rew.negative, rew.stronger1550, rew.ratio_outside);
fprintf('REW ratio by p_c4 threshold (ratio denominator requires finite W1550 > 1e-8):\n');
disp(rew.by_p_c4);
fprintf('QSO/component cells with any posterior REW ratio outside [1,2]:\n');
disp(rew.ratio_outside_qso_components);
fprintf('\np_c4 threshold report (denominator = active cells in each component):\n');
disp(thresholds);
if report.pass_numerical
    fprintf('PASS: no internal saved-output inconsistency was found.\n');
else
    fprintf('FAIL/WARN: inspect the non-zero counts above before using this run scientifically.\n');
end
end

function check_shapes(m, names, nq, nc)
for k = 1:numel(names)
    s = size(m, names{k});
    if startsWith(names{k}, 'posterior_') && endsWith(names{k}, '_draws')
        assert(numel(s) == 3 && s(1) == nq && s(3) == nc, '%s has wrong shape.', names{k});
    elseif ~ismember(names{k}, {'min_z_c4s', 'max_z_c4s'})
        assert(isequal(s, [nq nc]), '%s has wrong shape.', names{k});
    else
        assert(numel(m.(names{k})) == nq, '%s has wrong length.', names{k});
    end
end
end

function r = check_posterior_algebra(m, active, nq, nc)
p0 = m.p_no_c4; p1 = m.p_c4L1; p2 = m.p_c4;
lp0 = m.log_priors_no_c4; lpC = m.log_priors_c4;
ll0 = m.log_likelihoods_no_c4; ll1 = m.log_likelihoods_c4L1; ll2 = m.log_likelihoods_c4L2;
post0 = m.log_posteriors_no_c4; post1 = m.log_posteriors_c4L1; post2 = m.log_posteriors_c4L2;
invalid = 0;
for q = 1:nq
    for c = 1:nc
        if ~active(q,c), continue; end
        terms = [ll0(q,c)+lp0(q,c), ll1(q,c)+lpC(q,c), ll2(q,c)+lpC(q,c)];
        normv = logsumexp3(terms);
        % The processor saves log_posteriors_* as unnormalised log-joint
        % scores (log prior + log likelihood).  Only p_no_c4/p_c4L1/p_c4
        % are normalised across the three competing models.
        invalid = invalid + ~all(isfinite([terms normv post0(q,c) post1(q,c) post2(q,c)]));
        invalid = invalid + ~close_enough([post0(q,c) post1(q,c) post2(q,c)], terms);
        invalid = invalid + ~close_enough([p0(q,c) p1(q,c) p2(q,c)], exp(terms - normv));
    end
end
r = struct('invalid', invalid, 'message', sprintf('%d posterior/probability algebra mismatches.', invalid));
end

function r = check_map_and_intervals(m, active, nq, nc, ndraw, block_qsos)
% A matfile property does not support linear indexing (``m.x(:)``).  Load
% these two small vectors first, then make their orientation explicit.
zmin = m.min_z_c4s; zmin = zmin(:);
zmax = m.max_z_c4s; zmax = zmax(:);
zmap = m.map_z_c4L2; nmap = m.map_N_c4L2; smap = m.map_sigma_c4L2;
invalid = 0; checked = 0;
draw_names = {'posterior_z_c4L2_draws', 'posterior_logN_c4L2_draws', 'posterior_sigma_c4L2_draws', ...
    'posterior_REW_1548_draws', 'posterior_REW_1550_draws'};
q16_names = {'posterior_z_c4L2_16', 'posterior_logN_c4L2_16', 'posterior_sigma_c4L2_16', ...
    'posterior_REW_1548_16', 'posterior_REW_1550_16'};
q84_names = {'posterior_z_c4L2_84', 'posterior_logN_c4L2_84', 'posterior_sigma_c4L2_84', ...
    'posterior_REW_1548_84', 'posterior_REW_1550_84'};
num_blocks = ceil(nq / block_qsos);
fprintf('Checking MAP values and posterior intervals (%d QSO blocks) ...\n', num_blocks);
for block_number = 1:num_blocks
    q0 = (block_number - 1) * block_qsos + 1;
    qi = q0:min(q0+block_qsos-1,nq);
    d = cellfun(@(x) m.(x)(qi,:,:), draw_names, 'UniformOutput', false);
    lo = cellfun(@(x) m.(x)(qi,:), q16_names, 'UniformOutput', false);
    hi = cellfun(@(x) m.(x)(qi,:), q84_names, 'UniformOutput', false);
    active_block = active(qi,:);
    checked = checked + nnz(active_block);
    map_ok = isfinite(zmap(qi,:)) & isfinite(nmap(qi,:)) & isfinite(smap(qi,:)) & ...
        zmap(qi,:) >= zmin(qi) & zmap(qi,:) <= zmax(qi) & smap(qi,:) > 0;
    invalid = invalid + nnz(active_block & ~map_ok);
    for v = 1:numel(d)
        finite_draws = reshape(all(isfinite(d{v}), 2), [numel(qi) nc]);
        pct = prctile(d{v}, [16 84], 2);
        expected_lo = reshape(pct(:,1,:), [numel(qi) nc]);
        expected_hi = reshape(pct(:,2,:), [numel(qi) nc]);
        interval_ok = finite_draws & isfinite(lo{v}) & isfinite(hi{v}) & lo{v} <= hi{v} & ...
            close_elements(lo{v}, expected_lo) & close_elements(hi{v}, expected_hi);
        invalid = invalid + nnz(active_block & ~interval_ok);
    end
    if block_number == num_blocks || mod(block_number, 10) == 0
        fprintf('  posterior intervals: block %d/%d\n', block_number, num_blocks);
    end
end
r = struct('checked', checked, 'invalid', invalid, 'message', sprintf('%d active cells checked; %d MAP/draw/interval mismatches.', checked, invalid));
end

function r = check_rew_draws(m, p_c4, active, nq, nc, thresholds, block_qsos, make_plots)
thresholds = thresholds(:);
nt = numel(thresholds);
r = struct('finite',0,'negative',0,'stronger1550',0,'ratio_outside',0);
by_threshold = zeros(nt, 6); % finite pairs, valid ratios, W1550>W1548, ratio<1, ratio>2, outside
outside_qso_components = zeros(0, 4); % result QSO row, component, ratio<1 draws, ratio>2 draws
ratio_plot = [];
num_blocks = ceil(nq / block_qsos);
fprintf('Checking posterior REW draws (%d QSO blocks) ...\n', num_blocks);
for block_number = 1:num_blocks
    q0 = (block_number - 1) * block_qsos + 1;
    qi = q0:min(q0+block_qsos-1,nq);
    a = m.posterior_REW_1548_draws(qi,:,:); b = m.posterior_REW_1550_draws(qi,:,:);
    mask = repmat(reshape(active(qi,:), [numel(qi) 1 nc]), [1 size(a,2) 1]);
    x = a(mask); y = b(mask); good = isfinite(x) & isfinite(y);
    r.finite = r.finite + nnz(good); r.negative = r.negative + nnz(good & (x < -1e-10 | y < -1e-10));
    scale = max(1e-8, max(abs([x(good); y(good)]))); tol = 1e-7 * scale;
    r.stronger1550 = r.stronger1550 + nnz(good & y > x + tol);
    ratio = x(good & y > 1e-8) ./ y(good & y > 1e-8);
    r.ratio_outside = r.ratio_outside + nnz(ratio < 1-1e-5 | ratio > 2+1e-3);
    % Preserve the QSO/component identity of every affected spectrum.  For
    % this S_1_E_5000 result, result_qso_row is also all_quasar_ind used by
    % the processor's plotting path.
    for j = 1:numel(qi)
        for c = 1:nc
            if ~active(qi(j),c), continue; end
            xt = squeeze(a(j,:,c)); yt = squeeze(b(j,:,c));
            goodt = isfinite(xt) & isfinite(yt) & yt > 1e-8;
            ratiot = xt(goodt) ./ yt(goodt);
            low = nnz(ratiot < 1-1e-5); high = nnz(ratiot > 2+1e-3);
            if low + high > 0
                outside_qso_components(end+1,:) = [qi(j) c low high]; %#ok<AGROW>
            end
        end
    end
    for t = 1:nt
        selected_cells = active(qi,:) & p_c4(qi,:) >= thresholds(t);
        selected = repmat(reshape(selected_cells, [numel(qi) 1 nc]), [1 size(a,2) 1]);
        xt = a(selected); yt = b(selected); goodt = isfinite(xt) & isfinite(yt);
        ratiot = xt(goodt & yt > 1e-8) ./ yt(goodt & yt > 1e-8);
        finite_values = [xt(goodt); yt(goodt)];
        if isempty(finite_values), scale_t = 1e-8; else, scale_t = max(1e-8, max(abs(finite_values))); end
        tol_t = 1e-7 * scale_t;
        by_threshold(t,1) = by_threshold(t,1) + nnz(goodt);
        by_threshold(t,2) = by_threshold(t,2) + numel(ratiot);
        by_threshold(t,3) = by_threshold(t,3) + nnz(goodt & yt > xt + tol_t);
        by_threshold(t,4) = by_threshold(t,4) + nnz(ratiot < 1-1e-5);
        by_threshold(t,5) = by_threshold(t,5) + nnz(ratiot > 2+1e-3);
        by_threshold(t,6) = by_threshold(t,6) + nnz(ratiot < 1-1e-5 | ratiot > 2+1e-3);
    end
    if make_plots && numel(ratio_plot) < 2e6, ratio_plot = [ratio_plot; ratio(1:min(end,20000))]; end %#ok<AGROW>
    if block_number == num_blocks || mod(block_number, 10) == 0
        fprintf('  posterior REWs: block %d/%d\n', block_number, num_blocks);
    end
end
r.by_p_c4 = table(thresholds, by_threshold(:,1), by_threshold(:,2), by_threshold(:,3), ...
    by_threshold(:,4), by_threshold(:,5), by_threshold(:,6), ...
    by_threshold(:,3)./max(1,by_threshold(:,1)), by_threshold(:,6)./max(1,by_threshold(:,2)), ...
    'VariableNames', {'p_c4_threshold','finite_pairs','ratio_pairs','W1550_gt_W1548', ...
    'ratio_lt_1','ratio_gt_2','ratio_outside','W1550_gt_W1548_fraction','ratio_outside_fraction'});
r.ratio_outside_qso_components = array2table(outside_qso_components, 'VariableNames', ...
    {'result_qso_row','component','ratio_lt_1_draws','ratio_gt_2_draws'});
r.ratio_outside_qso_rows = unique(outside_qso_components(:,1), 'stable');
if make_plots && ~isempty(ratio_plot)
    figure('Color','w','Name','single-run C IV REW ratio'); histogram(ratio_plot, linspace(.5,2.5,101), 'Normalization','probability');
    xline(1,'k--'); xline(2,'k--'); xlim([.5 2.5]); grid on; xlabel('W_r(1548) / W_r(1550)'); ylabel('Fraction');
end
end

function r = check_sample_likelihoods(m, names, active, nq, nc, opt)
r = struct('checked',false,'mismatch',0,'message','Skipped by CheckSampleLikelihoods=false.');
if ~opt.CheckSampleLikelihoods, return; end
sample_name = 'sample_log_likelihoods_c4L2';
if ~ismember(sample_name,names), r.mismatch = inf; r.message = 'Saved L2 sample likelihoods are missing.'; return; end
stored = m.log_likelihoods_c4L2;
num_blocks = ceil(nq / opt.BlockQSOs);
fprintf('Checking L2 sample likelihoods (%d QSO blocks; streams the full sample cube) ...\n', num_blocks);
for block_number = 1:num_blocks
    q0 = (block_number - 1) * opt.BlockQSOs + 1;
    qi = q0:min(q0+opt.BlockQSOs-1,nq); block = m.(sample_name)(qi,:,:);
    for j = 1:numel(qi), for c = 1:nc
        if ~active(qi(j),c), continue; end
        x = squeeze(block(j,:,c));
        % Coverage-rejected proposals have log likelihood -Inf. Their zero
        % contributions must remain in the mean's proposal denominator.
        % NaN/+Inf in an active cell are invalid, rather than omitted samples.
        if any(isnan(x) | x == Inf)
            r.mismatch = r.mismatch + 1;
            continue;
        end
        mx = max(x);
        expected = mx + log(mean(exp(x-mx)));
        if ~isfinite(expected) || ~close_enough(stored(qi(j),c), expected), r.mismatch = r.mismatch + 1; end
    end, end
    if block_number == num_blocks || mod(block_number, 10) == 0
        fprintf('  L2 likelihoods: block %d/%d\n', block_number, num_blocks);
    end
end
r.checked = true; r.message = sprintf('%d L2 sample-to-marginal mismatches.', r.mismatch);
end

function r = plot_threshold_posteriors(m, p_c4, active, nq, nc, thresholds, block_qsos)
% Pool posterior draws by p_c4 threshold, retaining all saved draws.  The
% threshold masks select whole QSO/component cells, so the five quantities
% remain based on the same posterior samples.
draw_names = {'posterior_z_c4L2_draws', 'posterior_logN_c4L2_draws', ...
    'posterior_sigma_c4L2_draws', 'posterior_REW_1548_draws', 'posterior_REW_1550_draws'};
labels = {'z_{C IV}', 'log N_{C IV}', '\sigma_{C IV} (cm s^{-1})', ...
    'W_r(1548) (Angstrom)', 'W_r(1550) (Angstrom)'};
thresholds = thresholds(:)';
nt = numel(thresholds); nv = numel(draw_names);
num_blocks = ceil(nq / block_qsos);
lo = inf(1,nv); hi = -inf(1,nv); draw_counts = zeros(nt,nv);

% First pass: common finite range for all threshold curves of each quantity.
fprintf('Finding threshold-posterior histogram ranges (%d QSO blocks) ...\n', num_blocks);
for block_number = 1:num_blocks
    q0 = (block_number - 1) * block_qsos + 1;
    qi = q0:min(q0+block_qsos-1, nq);
    cell_mask = active(qi,:);
    d = cellfun(@(x) m.(x)(qi,:,:), draw_names, 'UniformOutput', false);
    for t = 1:nt
        selected_cells = cell_mask & p_c4(qi,:) >= thresholds(t);
        if ~any(selected_cells(:)), continue; end
        draw_mask = repmat(reshape(selected_cells, [numel(qi) 1 nc]), [1 size(d{1},2) 1]);
        for v = 1:nv
            x = d{v}(draw_mask); x = x(isfinite(x));
            if isempty(x), continue; end
            draw_counts(t,v) = draw_counts(t,v) + numel(x);
            % The least restrictive threshold determines shared edges.
            lo(v) = min(lo(v), min(x)); hi(v) = max(hi(v), max(x));
        end
    end
    if block_number == num_blocks || mod(block_number,10) == 0
        fprintf('  histogram ranges: block %d/%d\n', block_number, num_blocks);
    end
end

num_bins = 80; edges = cell(1,nv); counts = cell(nt,nv);
for v = 1:nv
    if ~isfinite(lo(v)) || ~isfinite(hi(v)), continue; end
    if lo(v) == hi(v), pad = max(1e-10, abs(lo(v))*1e-6); lo(v) = lo(v)-pad; hi(v) = hi(v)+pad; end
    edges{v} = linspace(lo(v), hi(v), num_bins+1);
    for t = 1:nt, counts{t,v} = zeros(1,num_bins); end
end

% Second pass: exact histogram counts, avoiding a large in-memory draw pool.
fprintf('Accumulating threshold-posterior histograms (%d QSO blocks) ...\n', num_blocks);
for block_number = 1:num_blocks
    q0 = (block_number - 1) * block_qsos + 1;
    qi = q0:min(q0+block_qsos-1, nq);
    cell_mask = active(qi,:);
    d = cellfun(@(x) m.(x)(qi,:,:), draw_names, 'UniformOutput', false);
    for t = 1:nt
        selected_cells = cell_mask & p_c4(qi,:) >= thresholds(t);
        if ~any(selected_cells(:)), continue; end
        draw_mask = repmat(reshape(selected_cells, [numel(qi) 1 nc]), [1 size(d{1},2) 1]);
        for v = 1:nv
            x = d{v}(draw_mask); x = x(isfinite(x));
            if ~isempty(x), counts{t,v} = counts{t,v} + histcounts(x, edges{v}); end
        end
    end
    if block_number == num_blocks || mod(block_number,10) == 0
        fprintf('  histogram counts: block %d/%d\n', block_number, num_blocks);
    end
end

colors = lines(nt);
draw_shape = size(m, 'posterior_z_c4L2_draws');
figure('Color','w','Name','C IV posterior samples by p_c4 threshold');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
legend_text = arrayfun(@(t) sprintf('p_{C IV} >= %.2f', t), thresholds, 'UniformOutput', false);
for v = 1:nv
    nexttile; hold on;
    centers = (edges{v}(1:end-1) + edges{v}(2:end))/2;
    bin_width = diff(edges{v});
    for t = 1:nt
        n = sum(counts{t,v});
        if n > 0, plot(centers, counts{t,v} ./ (n * bin_width), 'LineWidth', 1.3, 'Color', colors(t,:)); end
    end
    grid on; box on; xlabel(labels{v}); ylabel('Posterior-draw density');
    title(sprintf('%s: pooled posterior samples', labels{v}));
    if v == 1, legend(legend_text, 'Location','best'); end
end
nexttile; axis off;
text(0, .9, sprintf('All %d saved posterior draws per selected QSO/component.', draw_shape(2)), 'FontSize', 10);
for t = 1:nt
    text(0, .9-.12*t, sprintf('p_{C IV} >= %.2f: %d selected draws', thresholds(t), draw_counts(t,1)), 'Color', colors(t,:), 'FontSize', 10);
end
sgtitle('C IV posterior distributions by detection probability threshold');
r = struct('made',true,'thresholds',thresholds,'draw_counts',draw_counts, ...
    'bin_edges',{edges},'histogram_counts',{counts},'labels',{labels});
end

function t = threshold_table(p, active, thresholds)
thresholds = thresholds(:); nc = size(p,2); t = table(thresholds, 'VariableNames', {'p_c4_threshold'});
for c = 1:nc
    denom = nnz(active(:,c)); counts = arrayfun(@(x) nnz(active(:,c) & p(:,c) >= x), thresholds);
    t.(sprintf('component_%d_count',c)) = counts;
    t.(sprintf('component_%d_fraction',c)) = counts / max(1,denom);
end
end

function y = logsumexp3(x), mx = max(x); y = mx + log(sum(exp(x-mx))); end
function yes = close_enough(a,b), yes = all(abs(a(:)-b(:)) <= 1e-9 * max(1,max(abs([a(:);b(:)])))); end
function yes = close_elements(a,b), yes = abs(a-b) <= 1e-9 * max(1, max(abs(a), abs(b))); end
