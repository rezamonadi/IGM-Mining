% conditioned sampling of N and sigma_v given the W distribution
% using both PDF(W1) and PDF(W1)
% clc
clear
set_parameters;
% build_catalog;
load(sprintf('%s/catalog', processed_directory(release)), ...
    'EW1', 'EW2');
% load(sprintf('%s/catalog', processed_directory(release)), ...
%     'EW1');
EW1 = EW1*1e-8; % converting to cm from A
EW2 = EW2*1e-8; % converting to cm from A
min_EW1 = min(EW1);
max_EW1 = max(EW1);
min_EW2 = min(EW2);
max_EW2 = max(EW2);
PDF_EW1 = ksdensity(EW1, linspace(min_EW1, max_EW1, 1e3), 'Function', 'pdf');
PDF_EW2 = ksdensity(EW2, linspace(min_EW2, max_EW2, 1e3), 'Function', 'pdf');
max_PDF_EW1 = max(PDF_EW1);
max_PDF_EW2 = max(PDF_EW2);
min_PDF_EW1 = min(PDF_EW1);
min_PDF_EW2 = min(PDF_EW2);


m_e = 9.10938356e-28;
c   = 2.99792458e+10;
e   = 4.803204672997660e-10; 
fa = 0.189900;
fb  = 0.094750;
lambda_a = 1548e-8;
lambda_b = 1551e-8;
gamma_a = 2.643e+08;
gamma_b = 2.628e+08;
num_lambdas= 2000;             
lambda_array = linspace(1500e-8, 1600e-8, num_lambdas);                
d_lambda = lambda_array(2) - lambda_array(1);
nciv_samples = zeros(num_C4_samples,1);
sigma_samples = zeros(num_C4_samples,1);

% uniform_max_log_nciv = 17; 
% uniform_min_log_nciv= 12;
jj=0;
cc=0;
tic;
figure()
while jj<num_C4_samples
    cc=cc+1;
    % choosing random N and sigma_v
    N = rand*(uniform_max_log_nciv- uniform_min_log_nciv) + uniform_min_log_nciv;
    N = 10^N;
    sigma_v = rand*(max_sigma-min_sigma) +min_sigma;
    DLambda_a = sqrt(2)*sigma_v/c*lambda_a;
    DLambda_b = sqrt(2)*sigma_v/c*lambda_b;
    % finding phi_lambda and tau
    tau_a = zeros(num_lambdas,1);
    tau_b = zeros(num_lambdas,1);
    for k=1:num_lambdas
        % line a (1548A)
        I_phi_a = @(l) gamma_a*lambda_a^2/(4*pi*c)./((l-lambda_a).^2 +(gamma_a*lambda_a^2/4/pi/c)^2)*...
                   exp(-((l-lambda_array(k))/DLambda_a)^2)/sqrt(pi)/DLambda_a;
        phi_a = integral(I_phi_a, lambda_a-1000*DLambda_a, lambda_a+1000*DLambda_a, 'ArrayValued',true);
        tau_a(k,1) = N*e^2*lambda_a^2*fa/m_e/c^2*phi_a;

        % line b (1550A)
        I_phi_b = @(l) gamma_b*lambda_b^2/(4*pi*c)./((l-lambda_b).^2 +(gamma_b*lambda_b^2/4/pi/c)^2)*...
                   exp(-((l-lambda_array(k))/DLambda_b)^2)/sqrt(pi)/DLambda_b;
        phi_b = integral(I_phi_b, lambda_b-1000*DLambda_b, lambda_b+1000*DLambda_b, 'ArrayValued',true);
        tau_b(k,1) = N*e^2*lambda_b^2*fb/m_e/c^2*phi_b;
    end

    % calculating W_r having tau(N, sigma_v; lambda) -> Simpson integration
    EW_r_a=0;
    EW_r_b=0;
    
    % Simpson's Integration
    for k=1:num_lambdas-2
        % line 1548A
        f0 = 1-exp(-tau_a(k));
        fm = 1-exp(-tau_a(k+1));
        ff = 1-exp(-tau_a(k+2));
        EW_r_a = EW_r_a + d_lambda/3*(f0+ 4*fm + ff);

        % line 1550A
        f0 = 1-exp(-tau_b(k));
        fm = 1-exp(-tau_b(k+1));
        ff = 1-exp(-tau_b(k+2));
        EW_r_b = EW_r_b + d_lambda/3*(f0+ 4*fm + ff);
    end

    % evaluating PDF(EW1)    and PDF(EW2)
    % if (EW_r_a<=max_EW1 & EW_r_a>=min_EW1 & EW_r_b <=max_EW2 & EW_r_b>=min_EW2)
    % if (EW_r_a<=max_EW1 & EW_r_a>=min_EW1)
    Prob_EW_r_a = ksdensity(EW1, EW_r_a, 'Function', 'pdf');
    Prob_EW_r_b = ksdensity(EW2, EW_r_b, 'Function', 'pdf');
    if (Prob_EW_r_a>= (rand*(max_PDF_EW1 - min_PDF_EW1) + min_PDF_EW1) &...
        Prob_EW_r_b>= (rand*(max_PDF_EW2 - min_PDF_EW2) + min_PDF_EW2))
        jj=jj+1;
        fprintf('Iter: %5d, Sample %5d\n',cc, jj);
        nciv_samples(jj) = N;
        sigma_samples(jj) = sigma_v;
        EW_r_a_samples(jj) = EW_r_a;
        EW_r_b_samples(jj) = EW_r_b;
    end
end
toc;
rng('default');
sequence = scramble(haltonset(1), 'rr2');
% the first dimension can be used directly for the uniform prior over
% offsets
offset_z_samples  = sequence(1:num_C4_samples, 1)';
log_nciv_samples = log10(nciv_samples);
variables_to_save = {'offset_z_samples', 'sigma_samples',...
'log_nciv_samples', 'nciv_samples'};
save(sprintf('%s/civ_samples', processed_directory(training_release)), ...
     variables_to_save{:}, '-v7.3');

set(get(gca, 'XLabel'), 'String', '\sigma');
set(get(gca, 'YLabel'), 'String', 'N');
scatter(sigma_samples, log_nciv_samples, 10, 'r')
saveas(gcf, sprintf('Sigma-N-2L-2L-smp%d.png', num_C4_samples));
% clf 
figure()
ksdensity(EW1)
hold on
ksdensity(EW_r_a_samples)
hold on

% ksdensity(EW2)
% hold on
% ksdensity(EW_r_b_samples)

% legend('pdf(W1)','pdf(W1Sampled)')
legend('pdf(W1)','pdf(W1Sampled)') %, 'pdf(W2)', 'pdf(W2Sampled')
saveas(gcf, sprintf('PDF-Compare-smp%d-2L-nlambda%d.png', num_C4_samples, num_lambdas));