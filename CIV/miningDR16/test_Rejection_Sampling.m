% A sanity check for rejection sampling
clc
clear
num_xy_samples =1000;
max_x = 100; 
min_x= -100;
max_y = 100; 
min_y= -100;
jj=0;
cc=0;
tic;
r_dist = normrnd(15, 2, 1000,1);
max_r = max(r_dist);
min_r = min(r_dist);
max_PDF_r = max(ksdensity(r_dist));
min_PDF_r = min(ksdensity(r_dist));
x_samples = zeros(num_xy_samples,1);
y_samples = zeros(num_xy_samples,1);
r_samples = zeros(num_xy_samples,1);
while jj<num_xy_samples
    cc=cc+1;
    % choosing random N and sigma_v
    x_rand = rand*(max_x- min_x) + min_x;
    y_rand = rand*(max_y- min_y) + min_y;
    r_rand = sqrt(x_rand*x_rand + y_rand*y_rand);
    if (r_rand<=max_r & r_rand>=min_r)
        Prob_r = ksdensity(r_dist, r_rand, 'Function', 'pdf');
        if (Prob_r> (rand*(max_PDF_r - min_PDF_r) + min_PDF_r))
            jj=jj+1;
            fprintf('Iter: %15d, Sample %5d\n',cc, jj);
            x_samples(jj) = x_rand;
            y_samples(jj) = y_rand;
            r_samples(jj) = r_rand;
        end
    end
end
toc;
figure()
scatter(x_samples, y_samples, 10)
set(get(gca, 'XLabel'), 'String', 'x');
set(get(gca, 'YLabel'), 'String', 'y');

saveas(gcf, 'xy-test.png')
% clf 
figure()
ksdensity(r_dist)
hold on
ksdensity(r_samples)
% hold on

% ksdensity(EW2)
% hold on
% ksdensity(EW_r_b_samples)

legend('pdf(W1)','pdf(W1Sampled)')
% legend('pdf(W1)','pdf(W1Sampled)', 'pdf(W2)', 'pdf(W2Sampled')
saveas(gcf, 'PDF-Compare.png')