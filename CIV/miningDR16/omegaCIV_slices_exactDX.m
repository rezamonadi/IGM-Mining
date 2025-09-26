function [Omega_z, sigmaOmega_z, z_mid, out] = omegaCIV_slices_exactDX( ...
    phi_logN_z, var_phi_logN_z, logN_ctr, Delta_logN, ...
    z_edges, z_ranges, ...
    nSL, z_min_clipped, z_max_clipped, telluric_z, int_dX, opts)
% Ω_CIV in multiple redshift ranges using exact ΔX weighting per slice.
%
% Inputs
%   phi_logN_z      [nN x nZ]  completeness-corrected CDDF per dex per X
%   var_phi_logN_z  [nN x nZ]  variance of phi_logN_z (0 if unknown)
%   logN_ctr        [nN x 1]   log10 N centers (dex)
%   Delta_logN      scalar     width in dex (>0)
%   z_edges         [1 x (nZ+1)] z-bin edges used to build phi_logN_z
%   z_ranges        [nR x 2]   custom z-ranges [z_lo z_hi] (row per range)
%   nSL             scalar     number of sightlines
%   z_min_clipped   [nSL x 1]  per-sightline min usable z
%   z_max_clipped   [nSL x 1]  per-sightline max usable z
%   telluric_z      [nT x 2]   (can be empty []) masked z-intervals (applied to all SL)
%   int_dX          function   handle: dX = int_dX(zlo, zhi), same X as in phi
%   opts            struct     optional: opts.H0 (km/s/Mpc, default 70)
%
% Outputs
%   Omega_z         [nR x 1]   Ω_CIV per redshift range
%   sigmaOmega_z    [nR x 1]   1σ (diag-only) uncertainty from var_phi_logN_z
%   z_mid           [nR x 1]   midpoints of the requested ranges
%   out             struct     extras (K, H0_s, rho_c, W_range, X_range, phi_1D_multi, err_phi_1D_multi)

    % ---------- defaults & checks ----------
    if nargin < 12 || isempty(opts), opts = struct; end
    if ~isfield(opts,'H0'), opts.H0 = 70; end
    if isempty(var_phi_logN_z), var_phi_logN_z = zeros(size(phi_logN_z)); end
    assert(isscalar(Delta_logN) && Delta_logN > 0, 'Delta_logN must be a positive scalar.');

    [nN, nZ] = size(phi_logN_z);
    assert(all(size(var_phi_logN_z) == [nN, nZ]), 'var_phi_logN_z must match phi_logN_z size.');
    assert(numel(z_edges) == nZ+1, 'numel(z_edges) must be size(phi_logN_z,2)+1.');

    logN_ctr = logN_ctr(:);
    assert(numel(logN_ctr) == nN, 'logN_ctr length must equal size(phi_logN_z,1).');

    % scrub non-finite φ / Var(φ) safely
    phi_logN_z(~isfinite(phi_logN_z)) = 0;
    var_phi_logN_z(~isfinite(var_phi_logN_z) | var_phi_logN_z < 0) = 0;

    % ---------- constants (cgs) ----------
    c    = 2.99792458e10;                     % cm/s
    G    = 6.67430e-8;                        % cgs
    mp   = 1.67262192369e-24;                 % g
    mC   = 12*mp;                             % g
    H0_s = opts.H0 * 1e5 / 3.085677581e24;    % s^-1
    rho_c = 3*H0_s^2/(8*pi*G);                % g/cm^3
    K = (H0_s * mC) / (c * rho_c);            % ≈ 1.6e-22 for H0=70

    % ---------- grids ----------
    z_edges  = z_edges(:).';
    z_ranges = double(z_ranges);
    nR       = size(z_ranges,1);
    Nlin     = 10.^logN_ctr;

    % ---------- exact ΔX in each requested z-range, per z-bin ----------
    W_range = zeros(nR, nZ);
    hasTel  = ~isempty(telluric_z);
    if hasTel, telluric_z = double(telluric_z); end

    for r = 1:nR
        zlo_r = z_ranges(r,1);  zhi_r = z_ranges(r,2);
        for j = 1:nZ
            zb_lo = z_edges(j);  zb_hi = z_edges(j+1);
            lo = max(zb_lo, zlo_r);
            hi = min(zb_hi, zhi_r);
            if hi <= lo, continue; end

            dX_sum = 0;
            for s = 1:nSL
                zlo = max(z_min_clipped(s), lo);
                zhi = min(z_max_clipped(s), hi);
                if zhi <= zlo, continue; end

                dX_s = int_dX(zlo, zhi);
                if hasTel
                    for t = 1:size(telluric_z,1)
                        tlo = max(zlo, telluric_z(t,1));
                        thi = min(zhi, telluric_z(t,2));
                        if thi > tlo
                            dX_s = dX_s - int_dX(tlo, thi);
                        end
                    end
                end
                if dX_s < 0, dX_s = 0; end   % clamp tiny negatives
                dX_sum = dX_sum + dX_s;
            end
            W_range(r,j) = dX_sum;
        end
    end

    X_range = sum(W_range, 2);                 % total ΔX per requested range

    % ---------- path-weighted φ(logN) in each requested range ----------
    phi_1D_multi     = zeros(nN, nR);
    err_phi_1D_multi = zeros(nN, nR);

    for r = 1:nR
        if X_range(r) <= 0
            warning('No path in z-range [%.2f, %.2f].', z_ranges(r,1), z_ranges(r,2));
            continue;
        end
        wj = (W_range(r,:).');                                   % [nZ x 1]
        phi_r = (phi_logN_z * wj) / X_range(r);                  % [nN x 1]
        var_r = (var_phi_logN_z * (wj.^2)) / (X_range(r)^2);     % [nN x 1]
        phi_1D_multi(:,r)     = phi_r;
        err_phi_1D_multi(:,r) = sqrt(var_r);
    end

    % ---------- integrate Ω per range ----------
    Omega_z      = nan(nR,1);
    sigmaOmega_z = nan(nR,1);
    for r = 1:nR
        if X_range(r) <= 0, continue; end
        integ     = phi_1D_multi(:,r)     .* Nlin * Delta_logN;   % φ * N * dlog10 N
        integ_err = err_phi_1D_multi(:,r) .* Nlin * Delta_logN;   % diag-only
        Omega_z(r)      = K * nansum(integ);
        sigmaOmega_z(r) = K * sqrt(nansum(integ_err.^2));
    end

    % ---------- pack ----------
    z_mid = 0.5 * sum(z_ranges, 2);
    out = struct('K',K,'H0_s',H0_s,'rho_c',rho_c, ...
                 'W_range',W_range,'X_range',X_range, ...
                 'phi_1D_multi',phi_1D_multi,'err_phi_1D_multi',err_phi_1D_multi);
end
