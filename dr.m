clc
clear
uniform_max_log_nciv = 16.5;
set_parameters;
build_catalog;
% generate_c4_samples;
load(sprintf('%s/civ_samples-%s', processed_directory(training_release),...
            training_set_name), 'nciv_samples');

m_e = 9.10938356e-28;
c   = 2.99792458e+10;
e   = 4.803204672997660e-10; 
fa = 0.189900;
fb  = 0.094750;
lambda_a = 1548e-8;
lambda_b = 1551e-8;
gamma_a = 2.643e+08;
gamma_b = 2.628e+08;


b=10e5; 
% tau0 = 1.25393;
% % tau  = @(N, b, f, lambda) (0.7580*(N/1e13)*(f/0.4164)*(lambda/1215.7)*(10/b));
% tau = @(N, b, f, lambda) sqrt(pi)*e^2/m_e/c*N*f*lambda/b;

% W1 = @(N, b, f, lambda) sqrt(pi)*b/c*tau(N, b, f, lambda)/...
%                             (1+tau(N, b, f, lambda)/2/sqrt(2));
% W2 = @(N, b, f, lambda, gma) sqrt((2*b/c)^2*log(tau(N, b, f, lambda)/log(2))+...
%                    b/c*gma*lambda/c*(tau(N, b, f, lambda)-tau0)/sqrt(pi));

% DR = zeros(100,1);
% N_array = logspace(13.5,18, 100);
% for i=1:100
%     if (tau(N_array(i), b, fa, gamma_a)<=tau0)
%         DR(i) = W1(N_array(i), b, fa, lambda_a)/W1(N_array(i), b, fb, lambda_b);
%     else
%         DR(i) = W2(N_array(i), b, fa, lambda_a, gamma_a)/...
%                 W2(N_array(i), b, fb, lambda_b, gamma_b);
%     end
% end

% figure
% plot(N_array, DR)

sigma_v = b/sqrt(2);
% integrand of 
fa = @(nu, v) (1/sigma_v*exp(-v^2/2/sigma_v^2)*4*gamma_a/(16/pi^2/(nu - (1-v/c)*c/lambda_a)^2...
                                                                 + gamma_a^2) );

nu_array = linspace(c/lambda_a)                                                                 
for i=1:100                                                                 

