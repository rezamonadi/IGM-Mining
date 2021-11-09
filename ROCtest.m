clear
clc 
set_parameters
build_catalog
filename ='/home/reza/gpc/data/dr7/processed/processed_qsos_Rsigma-15-60.mat';

% filename = '/home/reza/gpc/data/dr7/processed/processed_qsos_RDLikelihood.mat';
load(filename);

ind_has_c4 =ismember(all_QSO_ID, c4_QSO_ID);
fig= figure();
y_score = p_c4;
y_true = ind_has_c4(test_ind);
[X,Y,T,AUC] =perfcurve(y_true, y_score, 'true');
plot(X,Y, 'lineWidth',5, 'Color', 'g')
legend(sprintf('AUC=%.5f',AUC))
set(get(gca, 'YLabel'), 'String', 'TPR');
set(get(gca, 'XLabel'), 'String', 'FPR');
set(gca, 'FontSize', 20)

testing_set_name = sprintf('sigma-%d-%d-N-%d-%d', min_sigma/1e5, max_sigma/1e5, uniform_min_log_nciv*10 , uniform_max_log_nciv*10)
exportgraphics(fig, sprintf('ROC-%s.pdf', testing_set_name)...
                ,'ContentType','vector')
