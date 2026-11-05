function report = diagnose_civ_pilot(result_file, proposal_file, varargin)
%DIAGNOSE_CIV_PILOT Pre-flight validation of a saved C IV processing block.
%
% Example (from the gpCIVIC directory):
%   report = diagnose_civ_pilot( ...
%       'processed_qsos_tst_DR16_S_1_E_1000.mat', ...
%       '/path/to/civ_samples_N-1200-1800-Sigma-5-115-Num-20000.mat', ...
%       'MinSigma', 5e5, 'MaxSigma', 115e5, 'MakePlots', true);
%
% The result file must be a v7.3 MAT-file.  proposal_file is optional, but
% supplying it enables exact MAP-to-proposal checks for N, sigma, and z.
% This script checks numerical integrity and internal consistency; it does
% NOT establish scientific completeness/purity on new DR16 sightlines.

if nargin < 1 || isempty(result_file)
    result_file = 'processed_qsos_tst_DR16_S_1_E_1000.mat';
end
if nargin < 2
    proposal_file = '';
end

% proposal_file is optional.  Permit the natural MATLAB call
% diagnose_civ_pilot(result_file, 'MakePlots', true) by recognizing that
% the second input is a name-value option rather than a file path.
option_names = {'MinSigma', 'MaxSigma', 'BlockQSOs', 'MakePlots', 'UseParallel'};
if (ischar(proposal_file) || (isstring(proposal_file) && isscalar(proposal_file))) && ...
        any(strcmpi(char(proposal_file), option_names))
    varargin = [{char(proposal_file)}, varargin];
    proposal_file = '';
end

p = inputParser;
p.addParameter('MinSigma', [], @(x) isempty(x) || (isscalar(x) && isfinite(x) && x > 0));
p.addParameter('MaxSigma', [], @(x) isempty(x) || (isscalar(x) && isfinite(x) && x > 0));
% 100 QSOs is about 112 MB for one 20000-by-7 double block.  It keeps the
% HDF5 read count low while remaining small relative to ordinary RAM.
p.addParameter('BlockQSOs', 100, @(x) isscalar(x) && x >= 1 && x == floor(x));
p.addParameter('MakePlots', false, @(x) islogical(x) && isscalar(x));
p.addParameter('UseParallel', false, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

assert(isfile(result_file), 'Result file not found: %s', result_file);
if ~isempty(proposal_file)
    assert(isfile(proposal_file), 'Proposal file not found: %s', proposal_file);
end

m = matfile(result_file);
names = who(m);
rew1548 = choose_name(names, {'REW_1548_dr16_samples', 'W_r_1548_samples'});
rew1550 = choose_name(names, {'REW_1550_dr16_samples', 'W_r_1550_samples'});
assert(~isempty(rew1548) && ~isempty(rew1550), ...
    'Could not find both 1548 and 1550 REW arrays in %s.', result_file);

sz48 = size(m, rew1548);
sz50 = size(m, rew1550);
assert(isequal(sz48, sz50) && numel(sz48) == 3, ...
    'The two REW arrays must have the same 3-D shape [QSO, sample, component].');
nq = sz48(1); ns = sz48(2); nc = sz48(3);

fprintf('\nC IV pilot pre-flight: %s\n', result_file);
fprintf('REW shape: %d QSOs x %d samples x %d components\n', nq, ns, nc);
if nq ~= 1000 || ns ~= 20000 || nc ~= 7
    warning('Expected pilot shape 1000 x 20000 x 7; inspecting the shape actually present.');
end

% A direct, physically motivated EW check.  Isolated C IV 1548 has the
% stronger oscillator strength, hence W1548 >= W1550.  The ratio approaches
% 2 when optically thin and 1 when saturated; a small tolerance absorbs
% quadrature and finite-window effects.
stats = struct('finite', zeros(1,nc), 'nan', zeros(1,nc), 'negative', zeros(1,nc), ...
    'zero', zeros(1,nc), 'stronger1550', zeros(1,nc), 'ratioOutside', zeros(1,nc));
ratio_for_plot = [];
num_blocks = ceil(nq / opt.BlockQSOs);
block_stats = zeros(num_blocks, nc, 6);
ratio_by_block = cell(num_blocks, 1);
if opt.UseParallel
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('local');
    end
    fprintf('Scanning REWs in parallel with %d workers (%d QSO blocks).\n', ...
        pool.NumWorkers, num_blocks);
    parfor block = 1:num_blocks
        [block_stats(block,:,:), ratio_by_block{block}] = ...
            scan_rew_block(result_file, rew1548, rew1550, block, opt.BlockQSOs, nq, nc, opt.MakePlots);
    end
else
    for block = 1:num_blocks
        [block_stats(block,:,:), ratio_by_block{block}] = ...
            scan_rew_block(result_file, rew1548, rew1550, block, opt.BlockQSOs, nq, nc, opt.MakePlots);
    end
end
stats.finite = reshape(sum(block_stats(:,:,1), 1), 1, nc);
stats.nan = reshape(sum(block_stats(:,:,2), 1), 1, nc);
stats.negative = reshape(sum(block_stats(:,:,3), 1), 1, nc);
stats.zero = reshape(sum(block_stats(:,:,4), 1), 1, nc);
stats.stronger1550 = reshape(sum(block_stats(:,:,5), 1), 1, nc);
stats.ratioOutside = reshape(sum(block_stats(:,:,6), 1), 1, nc);
if opt.MakePlots
    ratio_for_plot = vertcat(ratio_by_block{:});
    ratio_for_plot = ratio_for_plot(1:min(end, 2e6));
end

sample_total = nq * ns;
fprintf('\nREW integrity by component (each denominator is %d):\n', sample_total);
fprintf(' comp       finite        NaN/Inf   negative   W1550>W1548  ratio outside [1,2]\n');
for c = 1:nc
    fprintf('%5d %12d %14d %10d %13d %20d\n', c, stats.finite(c), ...
        stats.nan(c), stats.negative(c), stats.stronger1550(c), stats.ratioOutside(c));
end

% Later components are expected to be NaN for QSOs whose iterative search
% stopped earlier.  Each processed QSO contributes exactly ns REWs, so only
% a finite count that is not divisible by ns indicates an interrupted write.
incomplete_components = find(mod(stats.finite, ns) ~= 0);
processed_qsos_per_component = stats.finite / ns;
fprintf('Processed QSOs by component: %s\n', mat2str(processed_qsos_per_component));
if ~isempty(incomplete_components)
    warning('REW data end partway through a QSO for component(s): %s', ...
        mat2str(incomplete_components));
end

map_report = check_map_arrays(m, names, nq, nc);
probability_report = check_probabilities(m, names, nq, nc);
likelihood_report = check_marginal_likelihoods(m, names, nq, nc, opt.BlockQSOs);
proposal_report = struct('checked', false, 'mapRowsChecked', 0, 'mapMismatch', 0, ...
    'zMismatch', 0, 'nMismatch', 0, 'sigmaMismatch', 0, 'message', 'No proposal file supplied.');
if ~isempty(proposal_file)
    proposal_report = check_map_provenance(m, names, proposal_file, ns, nq, nc, opt);
end

fatal = sum(stats.negative) + sum(stats.stronger1550) + sum(stats.ratioOutside) + ...
    map_report.invalidFinite + probability_report.invalid + likelihood_report.mismatch;
if proposal_report.checked
    fatal = fatal + proposal_report.mapMismatch;
end

report = struct;
report.result_file = result_file;
report.proposal_file = proposal_file;
report.shape = sz48;
report.rew = stats;
report.map = map_report;
report.probabilities = probability_report;
report.marginal_likelihoods = likelihood_report;
report.provenance = proposal_report;
report.processed_qsos_per_component = processed_qsos_per_component;
report.pass_internal = (fatal == 0 && isempty(incomplete_components));

fprintf('\nMAP arrays: %s\n', map_report.message);
fprintf('Model probabilities: %s\n', probability_report.message);
fprintf('Marginal likelihoods: %s\n', likelihood_report.message);
fprintf('MAP provenance: %s\n', proposal_report.message);
if report.pass_internal
    fprintf('\nPASS: no internal numerical or proposal/MAP inconsistencies were found.\n');
else
    fprintf('\nFAIL/WARN: inspect the non-zero counts above before full-scale processing.\n');
end

if opt.MakePlots && ~isempty(ratio_for_plot)
    figure('Color', 'w', 'Name', 'C IV EW-ratio diagnostic');
    % Use fixed bins: a few tiny W1550 values otherwise produce enormous
    % ratios, causing MATLAB's automatic histogram range to collapse the
    % physically useful [1,2] interval into a single bar.
    edges = linspace(0.5, 2.5, 101);
    counts = histcounts(ratio_for_plot, edges);
    fraction = counts / numel(ratio_for_plot);
    centers = (edges(1:end-1) + edges(2:end)) / 2;
    bar(centers, fraction, 1, 'FaceColor', [0.30 0.60 0.80], 'EdgeColor', 'none');
    grid on; box on;
    xline(1, 'k--'); xline(2, 'k--'); xlim([0.5 2.5]);
    ylim([0, max(0.01, 1.1 * max(fraction))]);
    xlabel('W_r(1548) / W_r(1550)'); ylabel('Fraction');
    in_range = 100 * sum(counts) / numel(ratio_for_plot);
    title(sprintf('C IV EW ratio: %.1f%% of sampled finite ratios shown', in_range));
end
end

function [block_stats, ratio_for_plot] = scan_rew_block(result_file, rew1548, rew1550, block, block_qsos, nq, nc, make_plots)
% A worker-local MAT-file handle is required for parfor.  Each worker reads
% a disjoint QSO block, so no file is ever written or shared for mutation.
m = matfile(result_file);
q0 = (block - 1) * block_qsos + 1;
qi = q0:min(q0 + block_qsos - 1, nq);
a = m.(rew1548)(qi, :, :);
b = m.(rew1550)(qi, :, :);
block_stats = zeros(1, nc, 6);
ratio_for_plot = [];
for c = 1:nc
    x = a(:,:,c); y = b(:,:,c);
    finite = isfinite(x) & isfinite(y);
    block_stats(1,c,1) = nnz(finite);
    block_stats(1,c,2) = nnz(~finite);
    block_stats(1,c,3) = nnz(finite & (x < -1e-10 | y < -1e-10));
    block_stats(1,c,4) = nnz(finite & (x == 0 | y == 0));
    finite_values = [x(finite); y(finite)];
    if isempty(finite_values)
        scale = 1e-8;
    else
        scale = max(1e-8, max(abs(finite_values)));
    end
    tol = 1e-7 * scale;
    block_stats(1,c,5) = nnz(finite & y > x + tol);
    ratio = x(finite & y > 1e-8) ./ y(finite & y > 1e-8);
    block_stats(1,c,6) = nnz(ratio < 1-1e-5 | ratio > 2+1e-3);
    if make_plots && numel(ratio_for_plot) < 200000
        ratio_for_plot = [ratio_for_plot; ratio(1:min(end, 2000))]; %#ok<AGROW>
    end
end
end

function r = check_probabilities(m, names, nq, nc)
required = {'p_no_c4', 'p_c4L1', 'p_c4'};
r = struct('checked', false, 'invalid', 0, 'message', 'Posterior probability arrays not saved.');
if ~all(ismember(required, names)), return; end
p0 = m.p_no_c4; p1 = m.p_c4L1; p2 = m.p_c4;
if ~isequal(size(p0), [nq nc]) || ~isequal(size(p1), [nq nc]) || ~isequal(size(p2), [nq nc])
    r.invalid = inf; r.message = 'Posterior probability arrays have wrong shapes.'; return;
end
finite = isfinite(p0) & isfinite(p1) & isfinite(p2);
active = isfinite(p0) | isfinite(p1) | isfinite(p2);
r.invalid = nnz(active & ~finite) + nnz(finite & (p0 < -1e-12 | p1 < -1e-12 | p2 < -1e-12 | ...
    p0 > 1+1e-12 | p1 > 1+1e-12 | p2 > 1+1e-12 | abs(p0+p1+p2-1) > 1e-10));
r.checked = true;
r.message = sprintf('%d active cells; %d invalid or non-normalized.', nnz(finite), r.invalid);
end

function r = check_marginal_likelihoods(m, names, nq, nc, block_qsos)
sample_name = choose_name(names, {'sample_log_likelihoods_c4L2'});
stored_name = choose_name(names, {'log_likelihoods_c4L2'});
r = struct('checked', false, 'mismatch', 0, 'message', 'Doublet sample/stored likelihood arrays not both saved.');
if isempty(sample_name) || isempty(stored_name), return; end
if ~isequal(size(m, stored_name), [nq nc])
    r.mismatch = inf; r.message = 'log_likelihoods_c4L2 has the wrong shape.'; return;
end
% Read contiguous QSO blocks.  Reading one 1-by-20000 vector at a time
% causes thousands of costly HDF5 reads from a v7.3 MAT-file.
for q0 = 1:block_qsos:nq
    qi = q0:min(q0 + block_qsos - 1, nq);
    ll_block = m.(sample_name)(qi, :, :);
    stored_block = m.(stored_name)(qi, :);
    for q_local = 1:numel(qi)
        for c = 1:nc
        ll = squeeze(ll_block(q_local, :, c));
        if ~any(isfinite(ll)), continue; end
        mx = max(ll, [], 'omitnan');
        expected = mx + log(mean(exp(ll(isfinite(ll)) - mx)));
        if ~isfinite(stored_block(q_local,c)) || ~is_close(stored_block(q_local,c), expected)
            r.mismatch = r.mismatch + 1;
        end
        end
    end
end
r.checked = true;
r.message = sprintf('%d C IV doublet marginal-likelihood mismatches.', r.mismatch);
end

function r = check_map_arrays(m, names, nq, nc)
map_names = {'map_z_c4L2', 'map_N_c4L2', 'map_sigma_c4L2'};
missing = map_names(~ismember(map_names, names));
r = struct('missing', {missing}, 'invalidFinite', 0, 'message', '');
if ~isempty(missing)
    r.message = sprintf('Missing required doublet MAP array(s): %s.', strjoin(missing, ', '));
    r.invalidFinite = inf;
    return;
end
for k = 1:numel(map_names)
    assert(isequal(size(m, map_names{k}), [nq nc]), '%s has the wrong shape.', map_names{k});
end
z = m.map_z_c4L2; n = m.map_N_c4L2; s = m.map_sigma_c4L2;
active = isfinite(z) | isfinite(n) | isfinite(s);
consistent = isfinite(z) & isfinite(n) & isfinite(s);
r.invalidFinite = nnz(active & ~consistent) + nnz(consistent & (z <= -1 | s <= 0));
r.message = sprintf('%d active doublet MAP cells; %d invalid/incomplete.', nnz(consistent), r.invalidFinite);
end

function r = check_map_provenance(m, result_names, proposal_file, ns, nq, nc, opt)
pm = matfile(proposal_file); pnames = who(pm);
offset_z = load_vector(pm, pnames, {'offset_z_samples'});
logN = load_vector(pm, pnames, {'log_nciv_samples', 'log_N_samples', 'N_samples'});
sigma = load_vector(pm, pnames, {'sigma_civ_samples', 'sigma_samples'});
if isempty(sigma)
    off_sigma = load_vector(pm, pnames, {'offset_sigma_samples'});
    if ~isempty(off_sigma) && ~isempty(opt.MinSigma) && ~isempty(opt.MaxSigma)
        sigma = opt.MinSigma + (opt.MaxSigma - opt.MinSigma) * off_sigma;
    end
end
r = struct('checked', false, 'mapRowsChecked', 0, 'mapMismatch', 0, ...
    'zMismatch', 0, 'nMismatch', 0, 'sigmaMismatch', 0, 'message', '');
if numel(offset_z) ~= ns || numel(logN) ~= ns || numel(sigma) ~= ns
    r.message = ['Proposal vectors unavailable or wrong length. Need offset_z_samples, ', ...
        'log_nciv_samples (or equivalent), and sigma_civ_samples; for offset_sigma_samples ', ...
        'also provide MinSigma and MaxSigma.'];
    return;
end
ll_name = choose_name(result_names, {'sample_log_likelihoods_c4L2'});
if isempty(ll_name)
    r.message = 'No sample_log_likelihoods_c4L2 stored; cannot verify selected MAP index.';
    return;
end
minz = m.min_z_c4s; maxz = m.max_z_c4s;
if numel(minz) ~= nq || numel(maxz) ~= nq
    r.message = 'min_z_c4s/max_z_c4s missing or wrong length; cannot reconstruct z proposals.';
    return;
end
zmap = m.map_z_c4L2; nmap = m.map_N_c4L2; smap = m.map_sigma_c4L2;
% As above, use block reads rather than one HDF5 request per QSO/component.
for q0 = 1:opt.BlockQSOs:nq
    qi = q0:min(q0 + opt.BlockQSOs - 1, nq);
    ll_block = m.(ll_name)(qi, :, :);
    for q_local = 1:numel(qi)
        q = qi(q_local);
        zcand = minz(q) + (maxz(q) - minz(q)) * offset_z;
        for c = 1:nc
        ll = squeeze(ll_block(q_local, :, c));
        if ~any(isfinite(ll)), continue; end
        [~, ind] = max(ll, [], 'omitnan');
        r.mapRowsChecked = r.mapRowsChecked + 1;
        zbad = ~isfinite(zmap(q,c)) || ~is_close(zmap(q,c), zcand(ind));
        nbad = ~isfinite(nmap(q,c)) || ~is_close(nmap(q,c), logN(ind));
        sbad = ~isfinite(smap(q,c)) || ~is_close(smap(q,c), sigma(ind));
        r.zMismatch = r.zMismatch + zbad; r.nMismatch = r.nMismatch + nbad; r.sigmaMismatch = r.sigmaMismatch + sbad;
        end
    end
end
r.mapMismatch = r.zMismatch + r.nMismatch + r.sigmaMismatch;
r.checked = true;
r.message = sprintf('%d MAP selections checked: z=%d, N=%d, sigma=%d mismatches.', ...
    r.mapRowsChecked, r.zMismatch, r.nMismatch, r.sigmaMismatch);
end

function x = load_vector(m, names, candidates)
name = choose_name(names, candidates); x = [];
if ~isempty(name), x = m.(name); x = x(:); end
end

function name = choose_name(names, candidates)
name = ''; hit = candidates(ismember(candidates, names));
if ~isempty(hit), name = hit{1}; end
end

function yes = is_close(a, b)
yes = abs(a-b) <= 1e-10 * max([1, abs(a), abs(b)]);
end
