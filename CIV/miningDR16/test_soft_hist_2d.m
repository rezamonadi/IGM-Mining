function test_soft_hist_2d()
% Hand-checkable test for soft_hist_2d
% Grid: xedges=[0,1,2], yedges=[0,1,2]  → 2x2 bins
% 3 events, 4 posterior samples each, placed to make bin fractions obvious.

format short g

% ---- Inputs --------------------------------------------------------------
p = [1.0; 0.5; 0.4];     % inclusion probabilities for 3 events
S = 4;                   % samples per event

% Bin edges (2x2 grid)
xedges = [0, 1, 2];
yedges = [0, 1, 2];

% Posterior samples:
% Event 1 (p=1.0): all samples in bin (x1,y1) → (0.2, 0.2) repeated
samplesX(1,:) = [0.2 0.2 0.2 0.2];
samplesY(1,:) = [0.2 0.2 0.2 0.2];

% Event 2 (p=0.5): half in (x1,y2), half in (x2,y1)
%   (x1,y2) → (0.2, 1.2), (x2,y1) → (1.2, 0.2)
samplesX(2,:) = [0.2 0.2 1.2 1.2];
samplesY(2,:) = [1.2 1.2 0.2 0.2];

% Event 3 (p=0.4): one sample in each bin
%   (x1,y1),(x1,y2),(x2,y1),(x2,y2)
samplesX(3,:) = [0.2 0.2 1.2 1.2];
samplesY(3,:) = [0.2 1.2 0.2 1.2];

% ---- Run your function ---------------------------------------------------
out = ...
    soft_hist_2d(p, samplesX, samplesY, xedges, yedges);
H2D = out.H2D;
sigma2D = out.sigma2D;
Hx  = out.Hx;
sigmaHx = out.sigmaHx;
Hy  = out.Hy;
sigmaHy = out.sigmaHy;


% ---- Expected values (by hand) ------------------------------------------
% 2D expected counts:
%   Event 1: q = [1,0; 0,0]
%   Event 2: q = [0,0.25; 0.25,0]
%   Event 3: q = [0.1,0.1; 0.1,0.1]
H2D_exp = [1.1  0.35;
           0.35 0.10];

% Per-bin variances Var = sum_k q*(1-q)
Var2D_exp = [0.09   0.2775;
             0.2775 0.09];
sigma2D_exp = sqrt(Var2D_exp);

% 1D marginals (sum rows/cols of H2D)
Hx_exp = [1.45; 0.45];
Hy_exp = [1.45; 0.45];

% 1D covariance via event-level qx, qy
% x-marginal:
%   Event 1: qx=[1,0] → contrib diag(qx)-qx*qx' = 0
%   Event 2: qx=[0.25,0.25] → [[0.1875, -0.0625]; [-0.0625, 0.1875]]
%   Event 3: qx=[0.2,0.2]   → [[0.16,   -0.04  ]; [-0.04,    0.16  ]]
Covx_exp = [0.1875+0.16,  -0.0625-0.04;
            -0.0625-0.04,  0.1875+0.16];
sigmaHx_exp = sqrt(diag(Covx_exp));

% y-marginal is symmetric here
Covy_exp = Covx_exp;
sigmaHy_exp = sqrt(diag(Covy_exp));

% ---- Checks --------------------------------------------------------------
tol = 1e-12;

check('H2D',      H2D,      H2D_exp,      tol);
check('sigma2D',  sigma2D,  sigma2D_exp,  tol);
check('Hx',       Hx,       Hx_exp,       tol);
check('Hy',       Hy,       Hy_exp,       tol);
check('sigmaHx',  sigmaHx,  sigmaHx_exp,  tol);
check('sigmaHy',  sigmaHy,  sigmaHy_exp,  tol);

% Also sanity-check that totals match sum(p) since all samples lie in-range
W_exp = sum(p);                 % total expected count = 1.9
sumH2D = sum(H2D(:));
assert(abs(sumH2D - W_exp) < tol, 'Total expected count mismatch.');

% Pretty print
disp('--- Computed vs Expected ---');
disp('H2D (computed):');   disp(H2D);       disp('H2D (expected):');   disp(H2D_exp);
disp('sigma2D (computed):'); disp(sigma2D); disp('sigma2D (expected):'); disp(sigma2D_exp);
disp('Hx (computed):');    disp(Hx);        disp('Hx (expected):');    disp(Hx_exp);
disp('Hy (computed):');    disp(Hy);        disp('Hy (expected):');    disp(Hy_exp);
disp('sigmaHx (computed):'); disp(sigmaHx); disp('sigmaHx (expected):'); disp(sigmaHx_exp);
disp('sigmaHy (computed):'); disp(sigmaHy); disp('sigmaHy (expected):'); disp(sigmaHy_exp);

disp('✅ All checks passed.');

end

% ---------- helper: tolerant array compare ----------
function check(name, got, exp, tol)
    if ~isequal(size(got), size(exp))
        error('%s size mismatch: got %s, expected %s', name, mat2str(size(got)), mat2str(size(exp)));
    end
    d = max(abs(got(:) - exp(:)));
    if d > tol
        error('%s mismatch: max abs diff = %.3g (tol=%.1g)', name, d, tol);
    end
end
