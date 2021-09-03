% Author: Ming-Feng Ho (jibincat)
% Editted by: Reza Monadi (SunFlower)
% Rejection sampling:
%   Proc
%   1) Choose proposal q, exists c > 0 such that c q(X) >= f(X)
%   2) Sample X ~ q, sample Y ~ Unif(0, c * q(X)), given X
%   3) if Y <= f(X), then Z = X; Otherwise back to step 2.
%
% Output Z ~ f.
%
% Target f(X) is the EW kernel density PDF derived from DR7 catalog. 
%
% Our proposal dist q is uniformly sampling in (N, sigma), but it's
% not necessarily giving us a uniform dist in EW. So we need to
% first build the proposal q's PDF by uniformly sample in (N, sigma)
% and then use it to do rejection sampling.
% And it is crucial to find a constant value such that the proposal
% q * constant is always larger than target f(X).

%%%%%%%%%%%%%%%%%%%% Build proposal dist q %%%%%%%%%%%%%%%%%%%%
clear
rng('default');
set_parameters;

num_C4_samples = 5000;

% assume normal distributions for EWs
% EW1 = normrnd(0.8, 0.4, [1 10000]);
% EW2 = normrnd(1, 0.3, [1 10000]);
EW_cat=load(sprintf('%s/catalog', processed_directory(release)), ...
    'EW1', 'EW2');
EW1 = EW_cat.EW1;
EW2 = EW_cat.EW2;
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

uniform_max_log_nciv = 16; 
uniform_min_log_nciv= 13;
jj=0;
tic;
figure()
while jj<num_C4_samples
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
    % W_r= trapz(lambda_array, (1-exp(-tau)));
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

    jj=jj+1;
    if mod(jj,10)==0
        fprintf('Sample %5d\n',jj);
    end
    nciv_samples(jj) = N;
    sigma_samples(jj) = sigma_v;
    EW_r_a_samples(jj) = EW_r_a;
    EW_r_b_samples(jj) = EW_r_b;
end
toc;

% clf 
figure()
ksdensity(EW1, 'Kernel', 'box')
hold on
ksdensity(EW_r_a_samples, 'Kernel', 'box')
hold on

ksdensity(EW2, 'Kernel', 'box')
hold on
ksdensity(EW_r_b_samples, 'Kernel', 'box')

% legend('pdf(W1)','pdf(W1Sampled)')
legend('pdf(W1)','pdf(W1Sampled)', 'pdf(W2)', 'pdf(W2Sampled')
saveas(gcf, sprintf('PDF-Compare-smp%d-2L.png', num_C4_samples));

variables_to_save = {'EW_r_a_samples', 'EW_r_b_samples'};
save(sprintf('%s/plain_EW_sampling', processed_directory(training_release)), ...
     variables_to_save{:}, '-v7.3');


% %%%%%%%%%%%%%%%%%%%% Rejection sampling %%%%%%%%%%%%%%%%%%%%
% clear
% rng('default');
% set_parameters;


% proposal_q = load(sprintf('%s/plain_EW_sampling', processed_directory(training_release)));

% % assume normal distributions for EWs
% EW_cat=load(sprintf('%s/catalog', processed_directory(release)), ...
%     'EW1', 'EW2');
% EW1 = EW_cat.EW1;
% EW2 = EW_cat.EW2;
% EW1 = EW1*1e-8; % converting to cm from A
% EW2 = EW2*1e-8; % converting to cm from A
% min_EW1 = min(EW1);
% max_EW1 = max(EW1);
% min_EW2 = min(EW2);
% max_EW2 = max(EW2);
% PDF_EW1 = ksdensity(EW1, linspace(min_EW1, max_EW1, 1e3), 'Function', 'pdf', 'Kernel', 'box');
% PDF_EW2 = ksdensity(EW2, linspace(min_EW2, max_EW2, 1e3), 'Function', 'pdf', 'Kernel', 'box');
% max_PDF_EW1 = max(PDF_EW1);
% max_PDF_EW2 = max(PDF_EW2);
% min_PDF_EW1 = min(PDF_EW1);
% min_PDF_EW2 = min(PDF_EW2);


% % Check if the proposal dist (q) times a constant can be larger than the PDF we try to sample
% % Say constant to be constant = 3.5
% constant_a = 2.5;
% figure()
% [q_a, xi] = ksdensity(proposal_q.EW_r_a_samples, linspace(min_EW1, max_EW1, 1e3), 'Function', 'pdf', 'Kernel', 'box');

% plot(xi, q_a * constant_a);
% hold on
% plot(xi, PDF_EW1);
% % legend('Proposal q * c','Target PDF f');
% % saveas(gcf, sprintf('Proposal-target-dist.png'));
% hold on

% constant_b = 2.5;
% [q_b, xi] = ksdensity(proposal_q.EW_r_b_samples, linspace(min_EW2, max_EW2, 1e3), 'Function', 'pdf', 'Kernel', 'box');
% hold on
% plot(xi, q_b * constant_b);
% hold on
% plot(xi, PDF_EW2);
% legend(sprintf('Proposal q_a * %.2f', constant_a), 'Target PDF(EW1)',...
%        sprintf('Proposal q_b * %.2f', constant_b), 'Target PDF(EW2)');
% saveas(gcf, sprintf('Proposal-target-dist.png'));


% m_e = 9.10938356e-28;
% c   = 2.99792458e+10;
% e   = 4.803204672997660e-10; 
% fa = 0.189900;
% fb  = 0.094750;
% lambda_a = 1548e-8;
% lambda_b = 1551e-8;
% gamma_a = 2.643e+08;
% gamma_b = 2.628e+08;
% num_lambdas= 2000;             
% lambda_array = linspace(1500e-8, 1600e-8, num_lambdas);                
% d_lambda = lambda_array(2) - lambda_array(1);
% nciv_samples = zeros(num_C4_samples,1);
% sigma_samples = zeros(num_C4_samples,1);

% jj=0;
% cc=0;
% tic;
% figure()
% while jj<num_C4_samples
%     cc=cc+1;
%     % choosing random N and sigma_v
%     N = rand*(uniform_max_log_nciv- uniform_min_log_nciv) + uniform_min_log_nciv;
%     N = 10^N;
%     sigma_v = rand*(max_sigma-min_sigma) +min_sigma;
%     DLambda_a = sqrt(2)*sigma_v/c*lambda_a;
%     DLambda_b = sqrt(2)*sigma_v/c*lambda_b;
%     % finding phi_lambda and tau
%     tau_a = zeros(num_lambdas,1);
%     % tau_b = zeros(num_lambdas,1);
%     for k=1:num_lambdas
%         % line a (1548A)
%         I_phi_a = @(l) gamma_a*lambda_a^2/(4*pi*c)./((l-lambda_a).^2 +(gamma_a*lambda_a^2/4/pi/c)^2)*...
%                    exp(-((l-lambda_array(k))/DLambda_a)^2)/sqrt(pi)/DLambda_a;
%         phi_a = integral(I_phi_a, lambda_a-1000*DLambda_a, lambda_a+1000*DLambda_a, 'ArrayValued',true);
%         tau_a(k,1) = N*e^2*lambda_a^2*fa/m_e/c^2*phi_a;

%         % line b (1550A)
%         I_phi_b = @(l) gamma_b*lambda_b^2/(4*pi*c)./((l-lambda_b).^2 +(gamma_b*lambda_b^2/4/pi/c)^2)*...
%                    exp(-((l-lambda_array(k))/DLambda_b)^2)/sqrt(pi)/DLambda_b;
%         phi_b = integral(I_phi_b, lambda_b-1000*DLambda_b, lambda_b+1000*DLambda_b, 'ArrayValued',true);
%         tau_b(k,1) = N*e^2*lambda_b^2*fb/m_e/c^2*phi_b;
%     end

%     % calculating W_r having tau(N, sigma_v; lambda) -> Simpson integration
%     % W_r= trapz(lambda_array, (1-exp(-tau)));
%     EW_r_a=0;
%     EW_r_b=0;
    
%     % Simpson's Integration
%     for k=1:num_lambdas-2
%         % line 1548A
%         f0 = 1-exp(-tau_a(k));
%         fm = 1-exp(-tau_a(k+1));
%         ff = 1-exp(-tau_a(k+2));
%         EW_r_a = EW_r_a + d_lambda/3*(f0+ 4*fm + ff);

%         % line 1550A
%         f0 = 1-exp(-tau_b(k));
%         fm = 1-exp(-tau_b(k+1));
%         ff = 1-exp(-tau_b(k+2));
%         EW_r_b = EW_r_b + d_lambda/3*(f0+ 4*fm + ff);
%     end

    
%     % Rejection sampling:
%     %   Proc
%     %   1) Choose proposal q, exists c > 0 such that c * q(X) >= f(X)
%     %   2) Sample X ~ q, sample Y ~ Unif(0, c * q(X)), given X
%     %   3) if Y <= f(X), then Z = X; Otherwise back to step 2.
%     %
%     % Output Z ~ f.
%     Y_a = rand * constant_a * ksdensity(proposal_q.EW_r_a_samples, EW_r_a, 'Function', 'pdf', 'Kernel', 'box');
%     target_f_a = ksdensity(EW1, EW_r_a, 'Function', 'pdf', 'Kernel', 'box');

%     Y_b = rand * constant_b * ksdensity(proposal_q.EW_r_b_samples, EW_r_b, 'Function', 'pdf', 'Kernel', 'box');
%     target_f_b = ksdensity(EW2, EW_r_b, 'Function', 'pdf', 'Kernel', 'box');

%     if (Y_a <= target_f_a & Y_b <= target_f_b)
%         jj=jj+1;
%         fprintf('Iter: %5d, Sample %5d\n',cc, jj);
%         nciv_samples(jj) = N;
%         sigma_samples(jj) = sigma_v;
%         EW_r_a_samples(jj) = EW_r_a;
%         EW_r_b_samples(jj) = EW_r_b;
            
%     end
% end
% toc;
% rng('default');
% sequence = scramble(haltonset(1), 'rr2');
% % the first dimension can be used directly for the uniform prior over
% % offsets
% offset_z_samples  = sequence(1:num_C4_samples, 1)';
% log_nciv_samples = log10(nciv_samples);
% variables_to_save = {'offset_z_samples', 'sigma_samples',...
% 'log_nciv_samples', 'nciv_samples'};
% save(sprintf('%s/civ_samples_WR', processed_directory(training_release)), ...
%      variables_to_save{:}, '-v7.3');


% scatter(sigma_samples, log_nciv_samples);

% set(get(gca, 'XLabel'), 'String', '\sigma');
% set(get(gca, 'YLabel'), 'String', 'N');

% saveas(gcf, sprintf('test-proposal-Sigma-N-2L-2L-smp%d.png', num_C4_samples));
% % clf 
% figure()
% [f_target, xi] = ksdensity(EW1, linspace(min_EW1, max_EW1, 1e3), 'Function', 'pdf');
% plot(xi, f_target)
% hold on
% [q_proposal, xi] = ksdensity(EW_r_a_samples, linspace(min_EW1, max_EW1, 1e3), 'Function', 'pdf');
% plot(xi, q_proposal)
% hold on
% [f_target, xi] = ksdensity(EW2, linspace(min_EW2, max_EW2, 1e3), 'Function', 'pdf');
% plot(xi, f_target)
% hold on
% [q_proposal, xi] = ksdensity(EW_r_b_samples, linspace(min_EW2, max_EW2, 1e3), 'Function', 'pdf');
% plot(xi, q_proposal)

% legend('pdf(W1)','pdf(W1Sampled)', 'pdf(W2)', 'pdf(W2Sampled)' )
% saveas(gcf, sprintf('test-proposal-PDF-Compare-smp%d-2L.png', num_C4_samples));
