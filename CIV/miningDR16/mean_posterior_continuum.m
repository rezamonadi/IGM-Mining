function f_post = mean_posterior_continuum(y,mu,M,d)
  % Calculates Posterior Mean Continuum for plotting purposes
  y = y - mu;
  [n, k] = size(M);

  d_inv = 1 ./ d;
  D_inv_y = d_inv .* y;
  D_inv_M = d_inv .* M;

  B = M' * D_inv_M;
  B(1:(k + 1):end) = B(1:(k + 1):end) + 1;
  
  
  L = chol(B);
  C = L \ (L' \ D_inv_M');

  K_inv_y = D_inv_y - D_inv_M * (C * y);

  f_post = mu + M*(M'*K_inv_y);