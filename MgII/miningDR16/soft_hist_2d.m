function out = soft_hist_2d(p, events, xedges, yedges)
% SOFT_HIST_2D
% per-event posterior samples for (x,y) + inclusion probabilities p_k.
%
% INPUTS
%   p       : [N x 1] inclusion probabilities p_k in [0,1]
%   events  : 1xN struct array, each with fields:
%               - x : [S_k x 1] or [1 x S_k] posterior samples (or support pts) in x
%               - y : [S_k x 1] or [1 x S_k] posterior samples (or support pts) in y
%               - w : (optional) [S_k x 1] nonnegative weights for the samples
%   xedges  : [Bx+1 x 1] x bin edges
%   yedges  : [By+1 x 1] y bin edges
%
% OUTPUT (struct)
%   out.H2D        : [Bx x By] expected counts per 2D bin
%   out.sigma2D    : [Bx x By] per-bin std. dev. of realized counts
%   out.Hx_all     : [Bx x 1]  1D marginal over x (sum over y)
%   out.sigmaHx_all: [Bx x 1]  std. dev. for Hx_all (with covariance propagation)
%   out.Hy_all     : [By x 1]  1D marginal over y (sum over x)
%   out.sigmaHy_all: [By x 1]  std. dev. for Hy_all (with covariance propagation)
%   out.Hx_byY     : [Bx x By] "1D over x" for each y-bin (just H2D columns)
%   out.sigmaHx_byY: [Bx x By] per-bin std. dev. for those (equals sigma2D)
%   out.Hy_byX     : [By x Bx] "1D over y" for each x-bin (H2D rows, transposed)
%   out.sigmaHy_byX: [By x Bx] per-bin std. dev. for those (equals sigma2D')

    % ------------ input checks ------------
    p   = p(:);
    N   = numel(p);
    assert(isstruct(events) && numel(events)==N, 'events must be 1xN struct array matching p.');
    xedges = xedges(:);
    yedges = yedges(:);
    Bx = numel(xedges) - 1;
    By = numel(yedges) - 1;
    assert(Bx >= 1 && By >= 1, 'Need at least one bin along each axis.');

    % ------------ accumulators ------------
    H2D     = zeros(Bx, By);      % expected counts per 2D bin
    Var2D   = zeros(Bx, By);      % per-bin variances
    Hx_all  = zeros(Bx, 1);       % 1D marginal over x (sum over y)
    Hy_all  = zeros(By, 1);       % 1D marginal over y (sum over x)
    Covx    = zeros(Bx, Bx);      % covariance for Hx_all
    Covy    = zeros(By, By);      % covariance for Hy_all

    % ------------ main loop over events ------------
    for k = 1:N
        % pull samples / support
        xk = events(k).x(:);
        yk = events(k).y(:);
        assert(numel(xk) == numel(yk), 'events(k).x and events(k).y must have same length.');
        Sk = numel(xk);

        % per-sample weights (normalize to sum 1 over provided samples)
        if isfield(events(k), 'w') && ~isempty(events(k).w)
            wk = events(k).w(:);
            assert(numel(wk) == Sk, 'events(k).w must match number of samples.');
            wk = max(wk, 0);
            sw = sum(wk);
            if sw > 0
                wk = wk / sw;
            else
                wk = ones(Sk,1) / Sk;
            end
        else
            wk = ones(Sk,1) / Sk;
        end

        % per-event membership probabilities g_k(i,j)
        gk = membership2d(xk, yk, wk, xedges, yedges); % [Bx x By], sum <= 1

        % combine with inclusion probability
        qk = p(k) * gk; % [Bx x By]

        % accumulate expected counts and per-bin variances
        H2D   = H2D + qk;
        Var2D = Var2D + qk .* (1 - qk);

        % 1D projections with covariance propagation
        qx = sum(qk, 2);           % [Bx x 1], sum over y
        qy = sum(qk, 1).';         % [By x 1], sum over x
        Hx_all = Hx_all + qx;
        Hy_all = Hy_all + qy;
        Covx = Covx + diag(qx) - (qx * qx.');   % adds both var and cov terms
        Covy = Covy + diag(qy) - (qy * qy.');
    end

    % ------------ finalize ------------
    sigma2D     = sqrt(Var2D);
    sigmaHx_all = sqrt(diag(Covx));
    sigmaHy_all = sqrt(diag(Covy));

    % Pack also the per-slice "1D" views (these are just slices of H2D)
    out = struct();
    out.H2D          = H2D;
    out.sigma2D      = sigma2D;

    out.Hx_all       = Hx_all;
    out.sigmaHx_all  = sigmaHx_all;
    out.Hy_all       = Hy_all;
    out.sigmaHy_all  = sigmaHy_all;

    out.Hx_byY       = H2D;            % each column j is the 1D-over-x for y-bin j
    out.sigmaHx_byY  = sigma2D;
    out.Hy_byX       = H2D.';          % each column i is the 1D-over-y for x-bin i
    out.sigmaHy_byX  = sigma2D.';
end

% ---------- helper: per-event membership probabilities ----------
function g = membership2d(x, y, w, xedges, yedges)
% MEMBERSHIP2D  Return per-bin membership probabilities for a single event.
% Works on all MATLAB versions (no reliance on 'Weights' for histcounts2).
% Uses discretize + accumarray to form a weighted 2D histogram.
%
% Inputs:
%   x, y   : column vectors of samples or support points (same length)
%   w      : nonnegative weights (will be renormalized to sum 1 inside caller)
%   xedges : x bin edges (length Bx+1)
%   yedges : y bin edges (length By+1)
%
% Output:
%   g      : [Bx x By] matrix with entries in [0,1], sum(g(:)) <= 1

    Bx = numel(xedges) - 1;
    By = numel(yedges) - 1;

    % Bin indices (NaN if outside the range)
    ix = discretize(x, xedges);
    iy = discretize(y, yedges);

    % Keep only samples that fall inside the 2D grid
    m  = ~isnan(ix) & ~isnan(iy) & (w > 0);
    if ~any(m)
        g = zeros(Bx, By);
        return;
    end

    ix = ix(m);
    iy = iy(m);
    ww = w(m);

    % Accumulate weights into (ix, iy) cells
    g = accumarray([ix, iy], ww, [Bx, By], @sum, 0);

    % Note: weights were already normalized in the caller, so g is a
    % membership-probability table for this event (sums to <= 1).
end
