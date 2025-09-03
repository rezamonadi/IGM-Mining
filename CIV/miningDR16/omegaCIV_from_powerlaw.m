function [Omega, sigmaOmega, out] = omegaCIV_from_powerlaw(beta, sigma_beta, f0, sigma_f0, N0, Nmin, Nmax, opts)
% Compute Omega_CIV and its 1σ error for a power-law CDDF:
%   f(N) = f0 * (N/N0)^(-beta),  N in cm^-2, f in (number / cm^2 / dX).
%
% Inputs
%   beta, sigma_beta : slope and its 1σ uncertainty
%   f0, sigma_f0     : normalization at N0 and its 1σ uncertainty (set sigma_f0=0 if unknown)
%   N0               : pivot column density (cm^-2)
%   Nmin, Nmax       : integration limits (cm^-2)
%   opts (optional)  : struct with fields
%       .corr (default 0): correlation between ln f0 and (-beta) from your fit
%       .H0   (default 70): H0 in km/s/Mpc
%
% Outputs
%   Omega        : dimensionless Omega_CIV
%   sigmaOmega   : 1σ uncertainty propagated from beta (and f0 if provided)
%   out          : struct with pieces (K factor, integral, derivatives, etc.)
%
% Notes:
% - Error propagation is done in log-space:
%     ln Omega = ln K + ln f0 + beta*ln N0 + ln I(beta)
%   so Var[ln Omega] = Var[ln f0] + (d/dβ ln Omega)^2 Var[β] + 2*corr*σ_ln f0*σ_β*(-1)
%   where d/dβ ln Omega = ln N0 + d/dβ ln I(beta).
% - If you don’t have sigma_f0, pass 0 (only beta error is used).

    if nargin < 8 || isempty(opts), opts = struct; end
    if ~isfield(opts,'corr'), opts.corr = 0; end
    if ~isfield(opts,'H0'),   opts.H0   = 70; end

    % ----- constants (cgs) -----
    G  = 6.67430e-8;                  % cm^3 g^-1 s^-2
    c  = 2.99792458e10;               % cm s^-1
    mp = 1.67262192369e-24;           % g
    mC = 12.0 * mp;                   % g  (ion mass ~ carbon nucleus)
    H0 = opts.H0 * 3.240779289e-20;   % s^-1  (km/s/Mpc -> s^-1)

    rho_c = 3*H0^2 / (8*pi*G);        % g cm^-3
    K = (H0*mC) / (c*rho_c);          % dimensionless prefactor

    % ----- integral over CDDF -----
    [I, dlogI_dBeta] = I_and_dlog(beta, Nmin, Nmax);

    % ----- Omega -----
    Omega = K * f0 * (N0^beta) * I;

    % ----- error propagation (log-space) -----
    if isempty(sigma_f0), sigma_f0 = 0; end
    sigma_ln_f0 = (sigma_f0>0) * (sigma_f0 / f0); % 0 if sigma_f0==0
    dlnOmega_dBeta = log(N0) + dlogI_dBeta;       % ∂ lnΩ / ∂β

    % Covariance between ln f0 and β: cov = corr * σ_ln_f0 * σ_β * (-1) [since fit is for -β]
    cov_lnfo_beta = -opts.corr * sigma_ln_f0 * sigma_beta;

    var_lnOmega = (sigma_ln_f0)^2 + (dlnOmega_dBeta*sigma_beta)^2 + 2*cov_lnfo_beta*dlnOmega_dBeta;
    var_lnOmega = max(var_lnOmega, 0);            % guard against tiny negatives
    sigmaOmega  = Omega * sqrt(var_lnOmega);

    % ----- extras -----
    out = struct('K',K,'I',I,'dlogI_dBeta',dlogI_dBeta,'rho_c',rho_c,'H0_sinv',H0);
end

function [I, dlogI_dBeta] = I_and_dlog(beta, Nmin, Nmax)
% I(beta) = ∫_{Nmin}^{Nmax} N^{1-β} dN
% dlogI_dBeta = d/dβ [ln I(beta)]
    s = 2 - beta;
    % Integral
    if abs(s) > 1e-8
        I = (Nmax.^s - Nmin.^s) ./ s;
        % derivative via analytic dI/dβ = -dI/ds
        A = Nmax.^s; B = Nmin.^s; lnA = log(Nmax); lnB = log(Nmin);
        dIds = ((A.*lnA - B.*lnB).*s - (A - B)) ./ (s.^2);
        dIdb = -dIds;
        dlogI_dBeta = dIdb ./ I;
    else
        % Near β=2, use stable limits
        I = log(Nmax./Nmin);  % limit as s->0
        % Numerical derivative for safety
        epsb = 1e-6;
        Iplus  = integral_N1mbeta(beta+epsb, Nmin, Nmax);
        Iminus = integral_N1mbeta(beta-epsb, Nmin, Nmax);
        dlogI_dBeta = (log(Iplus) - log(Iminus)) / (2*epsb);
    end
end

function I = integral_N1mbeta(beta, Nmin, Nmax)
    s = 2 - beta;
    if abs(s) > 1e-12
        I = (Nmax.^s - Nmin.^s) ./ s;
    else
        I = log(Nmax./Nmin);
    end
end
