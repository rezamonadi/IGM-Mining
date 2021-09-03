% Confusion matrix builder 

clear
set_parameters;
build_catalog;
training_set_name = 'Rvoigt-10000Smp0-tr-80-20-vCutKathy-3000'
% training_set_name = 'Rvoigt-10000Smp0-tr-80-20-vCutReza-5000'

filename = sprintf('%s/processed_qsos_R%s', ...
    processed_directory(release), ...
    training_set_name);
load(filename, 'test_ind', 'p_L1','p_no_c4', 'map_z_c4L2', 'training_set_name',...
                'map_N_c4L2');
cooksey_catalog = ...
fitsread('data/dr7/distfiles/dr7qso_CIV_noBAL.fit', 'binarytable');
all_zqso          = cooksey_catalog{4};
ind_has_c4 =ismember(all_QSO_ID, c4_QSO_ID);

a = 1e-5;
p_c4 = 1- p_no_c4 -a*p_L1;
tr=0.85;

TP = nnz(ind_has_c4(test_ind) & p_c4>tr);
TN = nnz(~ind_has_c4(test_ind) & p_c4<tr);
FN = nnz(ind_has_c4(test_ind) & p_c4<tr );
ind_FP = ~ind_has_c4(test_ind) & p_c4>tr;
FP = nnz(ind_FP);
P = nnz(ind_has_c4(test_ind)  );
N = nnz(~ind_has_c4(test_ind));
confusion_matrix=[TP/P, FN/P; FP/N, TN/N];
Accuracy = (TP+TN)/(P+N);
ErrorRate = (FP+FN)/(P+N);
Sensitivity = TP/P;
Specificity = TN/N;
fprintf('p:%.2f\nFP:%d\nCM:[%.4f, %.4f\n    %.4f, %.4f]\nAccuracy:%.4f\nError Rate:%.4f\n'...
,tr, FP, TP/P, FN/P, FP/N,TN/N,Accuracy,ErrorRate);

% Checking if our FPs are intrinsic
% meaning that (Z_qso - MAP(Z_civ))/(1+Z_qso)>3000km/s/c

ind_intrinsic = (all_zqso(test_ind) - map_z_c4L2 )./(all_zqso(test_ind)+1)< kms_to_z(3000);
nnz(ind_intrinsic)