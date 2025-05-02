% % clear
% clc 
% cores= 4; 
% num_C4_samples = 20000;
% cataloging = 0;
% preloading = 0;
% learning   = 0;
% sampling   = 1;
% plotting   = 0;
% processing = 0;
% dv_mask    = 350;
% HPCC = 1;
% voigtPrep = 0;
% voigtType = 1;
% maskType = 1;
% priorType = 1;
% OccamRazor = 1;
% statTesting =0;
% MaskingProb = 0;
% set_parameters_dr7
% variables_to_load= {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7', ...
% 'all_QSO_ID', 'all_RA', 'all_DEC', 'all_zqso',...
% 'all_N_civ','all_z_civ3', 'all_RATING', 'c4_QSO_ID' };
% load(sprintf('%s/catalog', processed_directory(release)), ...
%     variables_to_load{:});
% % load('data/dr7/processed/processed_qsos_tst_OCCAML1-Sigma-4000000-nSamp-10000.mat');
% load('data/dr7/processed/processed_qsos_tst_1%-max-MgII-10.mat');
% load('data/dr7/processed/processed_qsos_tst_1%-norm2.mat');
load('data/dr7/processed/processed_qsos_tst_1%-norm4-20k.mat');


ind_has_MgII =ismember(all_QSO_ID, MgII_QSO_ID);
Z_Seyfeft = all_z_MgII3(test_ind,:);
z_qso = all_zqso(test_ind);
nTest = nnz(test_ind);
fig= figure();
% 
for quasar_ind=1:nTest
    DZ_min = 1000;
    i_min = -1;
% 
    p_MgII_best(quasar_ind)= p_MgII(quasar_ind,1);
    p_no_MgII_best(quasar_ind) = p_no_MgII(quasar_ind,1);
    p_MgIIL1_best(quasar_ind) = p_MgIIL1(quasar_ind,1);
    for i=1:10
        for j=1:17
            if (Z_Seyfeft(quasar_ind,j)>0)
                DZ = abs(Z_Seyfeft(quasar_ind, j) - map_z_MgIIL2(quasar_ind,i));
                if DZ<DZ_min
                    DZ_min = DZ;
                    i_min = i;
                end
            end
        end
    end
    if DZ_min< kms_to_z(dv_mask)*(1+z_qso(quasar_ind))
        p_MgII_best(quasar_ind)= p_MgII(quasar_ind, i_min);
        p_MgIIL1_best(quasar_ind) =p_MgIIL1(quasar_ind, i_min);
        p_no_MgII_best(quasar_ind) = p_no_MgII(quasar_ind,i_min);
% 
     end
end
fig = figure();
y_score = log(p_MgII_best./p_no_MgII_best);
y_true = ind_has_MgII(test_ind);
[X,Y,T,AUC] =perfcurve(y_true, y_score, 'true');
p=plot(X,Y, 'lineWidth',3);
p.Color = [0.3, 0.1, 0.5, 0.6];
legend(sprintf('AUC=%.3f',AUC),'location','southeast')
set(get(gca, 'YLabel'), 'String', 'True Positive Rate');
set(get(gca, 'XLabel'), 'String', 'False Positive Rate');
set(gca, 'FontSize', 14)
% 
% % testing_set_name = sprintf('sigma-%d-%d-N-%d-%d', min_sigma/1e5, max_sigma/1e5, uniform_min_log_nciv*10 , uniform_max_log_nciv*10)
exportgraphics(fig, 'ROC.png', 'Resolution', 800)
%                 % DZ = abs(Z_C13(quasar_ind, i) - map_z_c4L2(quasar_ind, j, 1));
%                 % if p_c4(quasar_ind, j)>=tr & DZ<kms_to_z(dv)*(1+zQSO_test(quasar_ind))
% 
% % fig = figure()
% % histogram(reshape(p_c4, length(p_c4)*7, 1), 40);
% % set(get(gca, 'YLabel'), 'String', 'Frequency');
% % set(get(gca, 'XLabel'), 'String', 'P(CIV)');
% % set(gca, 'FontSize', 15)
% % exportgraphics(fig, 'P_C4_hist.pdf', 'ContentType', 'vector')
