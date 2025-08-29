function out = soft_hist_2d(p, samplesX, samplesY, xedges, yedges)
% SOFTHIST2DROUTEA
% Expected 2D soft-count histogram and uncertainties from per-event posteriors.
%
% Inputs
%   p        : [N x 1] (or [1 x N]) inclusion probabilities p_k in [0,1]
%   samplesX : either [N x S] numeric array OR 1xN / Nx1 cell array; each cell/row has S samples of x for event k
%   samplesY : same layout as samplesX, S samples of y for each event k
%   xedges   : [Bx+1 x 1] vector of x bin edges (monotone ascending)
%   yedges   : [By+1 x 1] vector of y bin edges (monotone ascending)
%
% Outputs
%   H2D      : [Bx x By] expected counts per (x,y) bin
%   sigma2D  : [Bx x By] per-bin std. dev. for realized counts (incl. p_k + bin exclusivity)
%   Hx       : [Bx x 1]  x-marginal counts (summed over all y)
%   sigmaHx  : [Bx x 1]  std. dev. for x-marginal (uses proper bin-to-bin cov within each x row)
%   Hy       : [By x 1]  y-marginal counts (summed over all x)
%   sigmaHy  : [By x 1]  std. dev. for y-marginal (uses proper bin-to-bin cov within each y column)
%
% Math
%   For event k with posterior membership g_k(i,j) estimated from samples and inclusion p_k:
%     q_k(i,j) = p_k * g_k(i,j)
%   Then:
%     H2D(i,j)   = sum_k q_k(i,j)
%     Var2D(i,j) = sum_k q_k(i,j) * (1 - q_k(i,j))
%   x-marginal (vector over i):
%     qx_k(i) = sum_j q_k(i,j)
%     Cov_x  += diag(qx_k) - (qx_k * qx_k.')
%   y-marginal (vector over j):
%     qy_k(j) = sum_i q_k(i,j)
%     Cov_y  += diag(qy_k) - (qy_k * qy_k.')
%
% Notes
%   - Samples falling outside [xedges(1), xedges(end)) or [yedges(1), yedges(end)) are ignored
%     (i.e., treated as posterior mass outside the histogram window).
%   - If you want full 2D covariance, see the local helper at the end.

    % ---- validate & normalize inputs
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
            % No samples or zero inclusion prob → skip
            continue;
        end

        % 2D histogram of this event's posterior samples → counts per bin
        hk = histcounts2(xk(:), yk(:), xedges, yedges);   % [Bx x By]
        gk = hk / Sk;                                     % membership probabilities per (i,j)

        qk = p(k) * gk;                                   % per-bin success probs for realized count

        % 2D expected counts and per-bin variances
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
% Convert input samples to {1xN} cell arrays of row vectors for uniform handling.

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
        % Expect numeric arrays of size [N x S]
        assert(isnumeric(samplesX) && isnumeric(samplesY), 'samples must be numeric or cell.');
        assert(size(samplesX,1) == N && isequal(size(samplesX), size(samplesY)), ...
               'When numeric, samplesX and samplesY must be [N x S] with matching sizes.');
        S = size(samplesX,2);
        sxList = arrayfun(@(k) samplesX(k,1:S), 1:N, 'UniformOutput', false);
        syList = arrayfun(@(k) samplesY(k,1:S), 1:N, 'UniformOutput', false);
    end
end

