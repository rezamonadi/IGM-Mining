% Confusion matrix builder 

clear
set_parameters_dr7;
build_catalog_dr7;
filename = 'data/dr7/processed/processed_qsos_trn-mnp-200_tst-SnglMdl-1-f_L1_sigma-0.50-f_L1_N-0.10-Rsl-0-Rjc-0-nAvg-0-N-12.88d-15.80-sigma-10.00-50.00.mat';
variables_to_load = { 'training_set_name', ...
    'prior_ind', 'release', ...
    'test_ind', 'prior_z_qso_increase', ...
    'max_z_cut', 'min_z_c4s', 'max_z_c4s', ...
    'log_priors_no_c4', 'log_priors_c4', ...
    'log_likelihoods_no_c4',  ...
    'sample_log_likelihoods_c4L2', 'log_likelihoods_c4L2'...
    'log_posteriors_no_c4', 'log_posteriors_c4L1', 'log_posteriors_c4L2',...
    'model_posteriors', 'p_no_c4', ...
    'map_z_c4L2', 'map_N_c4L2', 'map_sigma_c4L2', 'p_c4_L2','p_c4_L1' 'EW1_spline', ...
    'EW2_spline', 'data_civ', 'EW1_fine', 'EW2_fine'};
load(filename, variables_to_load{:});
cooksey_catalog = ...
fitsread('data/dr7/distfiles/dr7qso_CIV_noBAL.fit', 'binarytable');
load(sprintf('%s/preloaded_qsos', processed_directory(release)));
all_zqso          = cooksey_catalog{4};
ind_has_c4 =ismember(all_QSO_ID_dr7, c4_QSO_ID);


tr=0.85;
predicted_c4s =  (p_c4_L2>tr);
ind_has_c4_tst = ind_has_c4(test_ind);
ind_TP = (ind_has_c4_tst & predicted_c4s);
TP = nnz(ind_TP );
TN = nnz(~ind_has_c4_tst & ~predicted_c4s);
FN = nnz(ind_has_c4_tst & ~predicted_c4s);
ind_FP = ~ind_has_c4_tst & predicted_c4s;
FP = nnz(ind_FP);
P = nnz(ind_has_c4_tst );
N = nnz(~ind_has_c4_tst);
confusion_matrix=[TP/P, FN/P; FP/N, TN/N];
Accuracy = (TP+TN)/(P+N);
ErrorRate = (FP+FN)/(P+N);
Sensitivity = TP/P;
Specificity = TN/N;
fprintf('tr:%.2f\nFP:%d\nCM:[%.3f, %.3f; %.3f, %.3f]\n'...
    ,tr, FP, TP/P, FN/P, FP/N,TN/N);%,Accuracy,ErrorRate);
fig= figure();
y_score = log(p_c4_L2./log(exp(p_no_c4)+exp(p_c4_L1)));
y_true = ind_has_c4(test_ind);
[X,Y,T,AUC1] =perfcurve(y_true, y_score, 'true');
plot(X,Y)
hold on
set(get(gca, 'YLabel'), 'String', 'TPR');
set(get(gca, 'XLabel'), 'String', 'FPR');
% set(get(gca, 'Title'), 'String', sprintf('p:%.2f, FP:%d\nCM:[%.4f, %.4f; %.4f, %.4f]\nAccuracy:%.4f, Error Rate:%.4f\n',...
% tr, FP, TP/P, FN/P, FP/N,TN/N,Accuracy,ErrorRate));

filename = 'data/dr7/processed/processed_qsos_trn-best-tr90_tst-SnglMdl-0-Rsl-0-Rjc-0-nAvg-0-N-1288-1580-sigma-10-50-nc4-50000-vCut-3000-prior-zqso-inc-10.mat';
variables_to_load = { 'training_set_name', ...
    'prior_ind', 'release', ...
    'test_ind', 'prior_z_qso_increase', ...
    'max_z_cut', 'min_z_c4s', 'max_z_c4s', ...
    'log_priors_no_c4', 'log_priors_c4', ...
    'log_likelihoods_no_c4',  ...
    'sample_log_likelihoods_c4L2', 'log_likelihoods_c4L2'...
    'log_posteriors_no_c4', 'log_posteriors_c4L1', 'log_posteriors_c4L2',...
    'model_posteriors', 'p_no_c4', ...
    'map_z_c4L2', 'map_N_c4L2', 'map_sigma_c4L2', 'p_c4','p_c4_L1' 'EW1_spline', ...
    'EW2_spline', 'data_civ', 'EW1_fine', 'EW2_fine'};
load(filename, variables_to_load{:});
cooksey_catalog = ...
fitsread('data/dr7/distfiles/dr7qso_CIV_noBAL.fit', 'binarytable');
load(sprintf('%s/preloaded_qsos', processed_directory(release)));
all_zqso          = cooksey_catalog{4};
ind_has_c4 =ismember(all_QSO_ID_dr7, c4_QSO_ID);

tr=0.85;
predicted_c4s =  (p_c4>tr);
ind_has_c4_tst = ind_has_c4(test_ind);
ind_TP = (ind_has_c4_tst & predicted_c4s);
TP = nnz(ind_TP );
TN = nnz(~ind_has_c4_tst & ~predicted_c4s);
FN = nnz(ind_has_c4_tst & ~predicted_c4s);
ind_FP = ~ind_has_c4_tst & predicted_c4s;
FP = nnz(ind_FP);
P = nnz(ind_has_c4_tst );
N = nnz(~ind_has_c4_tst);
confusion_matrix=[TP/P, FN/P; FP/N, TN/N];
Accuracy = (TP+TN)/(P+N);
ErrorRate = (FP+FN)/(P+N);
Sensitivity = TP/P;
Specificity = TN/N;
fprintf('tr:%.2f\nFP:%d\nCM:[%.3f, %.3f; %.3f, %.3f]\n'...
    ,tr, FP, TP/P, FN/P, FP/N,TN/N);%,Accuracy,ErrorRate);
y_score = log(p_c4./p_no_c4);
y_true = ind_has_c4(test_ind);
[X,Y,T,AUC2] =perfcurve(y_true, y_score, 'true');
plot(X,Y)
legend(sprintf('AUC=%.4f',AUC1), sprintf('AUC=%.4f', AUC2));


exportgraphics(fig, 'ROC-compare.pdf','ContentType','vector')
