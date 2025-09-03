function [beta,f0,beta_err,f0_err,stats] = fit_cddf_powerlaw(N,f,sigma_f,N0)
    N=N(:); f=f(:); s=sigma_f(:);
    ok = isfinite(N)&isfinite(f)&isfinite(s)&N>0&f>0&s>0; N=N(ok); f=f(ok); s=s(ok);
    if nargin<4||isempty(N0), w=(f./s).^2; N0=exp(sum(w.*log(N))/sum(w)); end
    x = log(N./N0); y = log(f); sy = s./f; w = 1./(sy.^2);
    mdl = fitlm(table(x,y,w),'y ~ 1 + x','Weights','w');
    a = mdl.Coefficients.Estimate(1);  b = mdl.Coefficients.Estimate(2);
    da = mdl.Coefficients.SE(1);       db = mdl.Coefficients.SE(2);
    f0 = exp(a);  beta = -b;  f0_err = f0*da;  beta_err = db;
    yh  = predict(mdl,table(x,y,w));
    chi2 = sum(((y - yh)./sy).^2);  dof = max(numel(y)-2,1);
    stats = struct('N0',N0,'chi2',chi2,'dof',dof,'chi2_red',chi2/dof, ...
                   'coef_cov', mdl.CoefficientCovariance);
end