% log_mvnpdf_low_rank: efficiently computes
%
%   log N(y; mu, MM' + diag(d))

function log_p = log_mvnpdf_low_rank(y, mu, M, d)

  log_2pi = 1.83787706640934534;

  [n, k] = size(M);

  y = y - mu;

  d_inv = 1 ./ d;
  D_inv_y = d_inv .* y;
  D_inv_M = d_inv .* M;
  % fprintf('S(D_inv_M)=%d-%d\n', size(D_inv_M));
  % use Woodbury identity, define
  %   B = (I + M' D^-1 M),
  % then
  %   K^-1 = D^-1 - D^-1 M B^-1 M' D^-1
  % fprintf('S(M)=%d-%d\n', size(M));
  % D_inv_M
  % M
  B = M' * D_inv_M;
  % B
  % fprintf('S(B)=%d-%d\n', size(B));
  % fprintf('isnan(M)=%d\n', nnz(isnan(M)));
  % fprintf('zeros(M)=%d\n', nnz((M==0)));
  % fprintf('inf(M)=%d\n', nnz((isinf(M))));

  % fprintf('isnan(d_inv)=%d\n', nnz(isnan(d_inv)));
  % fprintf('zeros(d_inv)=%d\n', nnz(d_inv==0))
  % fprintf('inf(d_inv)=%d\n', nnz(isinf(d_inv)))

  % fprintf('isnan(D_inv_M)=%d\n', nnz(isnan(D_inv_M)));
  % fprintf('zeros(D_inv_M)=%d\n', nnz(D_inv_M==0))
  % fprintf('inf(D_inv_M)=%d\n', nnz(isinf(D_inv_M)))

  % fprintf('isnan(B)=%d\n', nnz(isnan(B)));
  % fprintf('zeros(B)=%d\n', nnz(B==0))
  % fprintf('inf(B)=%d\n', nnz(isinf(B)))
  B(1:(k + 1):end) = B(1:(k + 1):end) + 1;
  
  
  L = chol(B);
  % C = B^-1 M' D^-1
  C = L \ (L' \ D_inv_M');

  K_inv_y = D_inv_y - D_inv_M * (C * y);

  log_det_K = sum(log(d)) + 2 * sum(log(diag(L)));

  log_p = -0.5 * (y' * K_inv_y + log_det_K + n * log_2pi);

end