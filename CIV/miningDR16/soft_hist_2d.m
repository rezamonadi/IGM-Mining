function out = soft_hist_2d(p, samplesX, samplesY, xedges, yedges, opts)
% SOFTHIST2DROUTEA
% Expected 2D soft-count histogram and uncertainties from per-event posteriors,
% with optional Poisson (shot) noise added to the variance.
%
% New option:
%   opts.addPoisson : logical (default=false). When true, adds Poisson variance:
%                     Var2D += H2D; Covx += diag(Hx); Covy += diag(Hy).

    % ---- validate & normalize inputs
    if nargin < 6 || isempty(opts), opts = struct; end
    if ~isfield(opts,'addPoisson'), opts.addPoisson = false; end

    p = p(:);
    N = numel(p);
    assert(N >= 1, 'p must have at least one event.');
    [sxList, syList] = normalizeSamples(samplesX, samplesY, N);

    Bx = numel(xedges) - 1;
    By = numel(yedges) - 1;
    assert(Bx >= 1 && By >= 1, 'Need at least one bin in x and y.');

    % ---- accumulators
    H2D   = zeros(Bx, By);
    Var2D = zeros(Bx, By);
    Covx  = zeros(Bx, Bx);
    Covy  = zeros(By, By);
    Hx    = zeros(Bx, 1);
    Hy    = zeros(By, 1);

    % ---- per-event contribution
    for k = 1:N
        xk = sxList{k};  yk = syList{k};
        assert(isvector(xk) && isvector(yk) && numel(xk) == numel(yk), ...
               'Event %d: x and y must be vectors with the same length.', k);
        Sk = numel(xk);
        if Sk == 0 || p(k) == 0
            continue;
        end

        % 2D histogram of this event's posterior samples → counts per bin
        hk = histcounts2(xk(:), yk(:), xedges, yedges);   % [Bx x By]
        gk = hk / Sk;                                     % membership probabilities per (i,j)

        qk = p(k) * gk;                                   % per-bin success probs for realized count

        % 2D expected counts and per-bin variances (bin-hopping/measurement term)
        H2D   = H2D   + qk;
        Var2D = Var2D + qk .* (1 - qk);

        % x-marginal contribution (vector over i), with covariance
        qxk = sum(qk, 2);                                 % [Bx x 1]
        Hx  = Hx + qxk;
        Covx = Covx + diag(qxk) - (qxk * qxk.');

        % y-marginal contribution (vector over j), with covariance
        qyk = sum(qk, 1).';                               % [By x 1]
        Hy  = Hy + qyk;
        Covy = Covy + diag(qyk) - (qyk * qyk.');
    end

    % ---- optional: add Poisson counting variance (AFTER accumulation)
    if opts.addPoisson
        % For unit weights, Poisson variance per 2D bin is Var = H2D.
        % For marginals, Poisson contributes only to the diagonal.
        Var2D = Var2D + H2D;          % per 2D bin
        Covx  = Covx  + diag(Hx);     % x-marginal bins (no cross terms)
        Covy  = Covy  + diag(Hy);     % y-marginal bins (no cross terms)
    end

    % ---- convert variances → standard deviations
    sigma2D = sqrt(Var2D);
    sigmaHx = sqrt(diag(Covx));
    sigmaHy = sqrt(diag(Covy));

    out.H2D      = H2D;
    out.sigma2D  = sigma2D;
    out.Hx       = Hx;
    out.sigmaHx  = sigmaHx;
    out.Hy       = Hy;
    out.sigmaHy  = sigmaHy;

    dx = diff(xedges);
    dy = diff(yedges);
    out.dx = dx;
    out.dy = dy;
    out.xcenters = xedges(1:end-1) + 0.5*dx;
    out.ycenters = yedges(1:end-1) + 0.5*dy;
end

% ---------- helpers ----------
function [sxList, syList] = normalizeSamples(samplesX, samplesY, N)
    if iscell(samplesX)
        assert(iscell(samplesY), 'samplesY must be a cell array if samplesX is.');
        assert(numel(samplesX) == N && numel(samplesY) == N, ...
               'Cell arrays must have N elements.');
        sxList = cell(1,N); syList = cell(1,N);
        for k = 1:N
            sx = samplesX{k}; sy = samplesY{k};
            sxList{k} = sx(:).';   % row
            syList{k} = sy(:).';
        end
    else
        assert(isnumeric(samplesX) && isnumeric(samplesY), 'samples must be numeric or cell.');
        assert(size(samplesX,1) == N && isequal(size(samplesX), size(samplesY)), ...
               'When numeric, samplesX and samplesY must be [N x S] with matching sizes.');
        S = size(samplesX,2);
        sxList = arrayfun(@(k) samplesX(k,1:S), 1:N, 'UniformOutput', false);
        syList = arrayfun(@(k) samplesY(k,1:S), 1:N, 'UniformOutput', false);
    end
end
