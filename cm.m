% Confusion matrix builder 


clear
set_parameters;
build_catalog;

% training_set_name ='-WNciv-to-1600';
filename = sprintf('%s/processed_qsos_R%s', ...
    processed_directory(release), ...
    training_set_name);
load(filename, 'test_ind', 'p_L1','p_no_c4', 'map_z_c4L2', 'training_set_name',...
                'map_N_c4L2');
% training_set_name
cooksey_catalog = ...
fitsread('data/dr7/distfiles/dr7qso_CIV_noBAL.fit', 'binarytable');
all_zqso          = cooksey_catalog{4};
ind_has_c4 =ismember(all_QSO_ID, c4_QSO_ID) & (filter_flags==0) & (all_RATING==3);

a = 1e-5;
p_c4 = 1- p_no_c4 -a*p_L1;
% p_c4 = p_c4(all_zqso(test_ind)> (map_z_c4L2+0.1));
% test_ind = test_ind(all_zqso(test_ind)> (map_z_c4L2+0.01));
min_N = 13.25;
% tr = 0.85
% for min_N=[13, 13.25, 13.5, 13.75, 14]
tr=0.9;
Dz = max_z_cut;
Dz = -100

TP = nnz(ind_has_c4(test_ind) & p_c4>tr  & (all_zqso(test_ind)> (map_z_c4L2+Dz)) );
TN = nnz(~ind_has_c4(test_ind) & p_c4<tr & (all_zqso(test_ind)> (map_z_c4L2+Dz)));
FN = nnz(ind_has_c4(test_ind) & p_c4<tr & (all_zqso(test_ind)> (map_z_c4L2+Dz)));
FP = nnz(~ind_has_c4(test_ind) & p_c4>tr & (all_zqso(test_ind)> (map_z_c4L2+Dz)) );
P = nnz(ind_has_c4(test_ind) & (all_zqso(test_ind)> (map_z_c4L2+Dz)) );
N = nnz(~ind_has_c4(test_ind)  & (all_zqso(test_ind)> (map_z_c4L2+Dz)));
confusion_matrix=[TP/P, FN/P; FP/N, TN/N];
Accuracy = (TP+TN)/(P+N);
ErrorRate = (FP+FN)/(P+N);
Sensitivity = TP/P;
Specificity = TN/N;
fprintf('p:%.2f\nFP:%d\nCM:[%.4f, %.4f\n    %.4f, %.4f]\nAccuracy:%.4f\nError Rate:%.4f\n'...
,tr, FP, TP/P, FN/P, FP/N,TN/N,Accuracy,ErrorRate);
% fprintf('min(N):%.2f\nFP:%d\n',min_N, FP);
fprintf('--------------------\n\n')


% tr=0.9;
% figure();
% % histogram(p_no_c4(ind_has_c4(test_ind) & p_c4<tr))
% % hold on
% histogram(p_L1(ind_has_c4(test_ind) & p_c4<0.1))
% % hold on
% % histogram(p_c4(ind_has_c4(test_ind) & p_c4<tr))
% % legend('pNull','pL1', 'pC4')

% xlabel('p(single line)')
% title('C13 detected a doublet but p(doublet)<0.1')
% saveas(gcf,'p_single.png')
