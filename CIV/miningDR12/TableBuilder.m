fprintf('loading...\n')
filename = sprintf('%s/processed_qsos_dr12_N-1250-1610-S-35-115-nc-10k.mat', processed_directory(releaseTest));
% load(filename);
load('CredIntervals/CIs95-10.mat');
load('REW_1548_DR12.mat');
load('ShortProcessedDR12.mat');
fid = fopen('Table95All_Pciv_positive_10.dat', 'w');
% test_ind = savingCat.test_ind;
% z_qsos =  zqso_dr12(test_ind);
% num_quasars    =     numel(z_qsos);
% QSO_ID_dr12 = QSO_ID_dr12(test_ind);
% map_z_c4L2 = savingCat.map_z_c4L2;
% map_N_c4L2 = savingCat.map_N_c4L2;
% map_sigma_c4L2 = savingCat.map_sigma_c4L2;
% p_c4 = savingCat.p_c4;
% p_c4L1 = savingCat.p_c4L1;

fprintf('Printing...\n')
for quasar_ind=1:10
    for num_c4=1:7
        % if (p_c4(quasar_ind, num_c4)>=0)
        if (p_c4(quasar_ind, num_c4)>=0)
            fprintf(fid, ...
                '%-18s  %.4f %.5f %.5f %.2f %.2f %03.2f %3.2f %.3f %.3f %.3f %.3f %.2f %.2f\n', ...
                all_QSO_ID_dr12{quasar_ind},...
                z_qsos(quasar_ind) ,...
                all_map_z_c4L2(quasar_ind, num_c4),...
                (CI_Z(quasar_ind,num_c4,2) - CI_Z(quasar_ind,num_c4,1))*0.5,...
                all_map_N_c4L2(quasar_ind, num_c4),...
                (CI_N(quasar_ind,num_c4,2)- CI_N(quasar_ind,num_c4,1))*0.5,...
                all_map_sigma_c4L2(quasar_ind, num_c4)/1e5,...
                (CI_Sigma(quasar_ind,num_c4,2)-CI_Sigma(quasar_ind,num_c4,1))*0.5/1e5,...
                REW_1548_DR12_voigt(quasar_ind, num_c4), ...
                (CI_W_r_1548(quasar_ind, num_c4,2)-CI_W_r_1548(quasar_ind, num_c4,1))*0.5, ...
                REW_1550_DR12_voigt(quasar_ind, num_c4), ...
                (CI_W_r_1550(quasar_ind, num_c4,2)-CI_W_r_1550(quasar_ind, num_c4,1))*0.5,...
                all_p_c4(quasar_ind, num_c4),...
                all_p_c4L1(quasar_ind, num_c4));
        end
    end
    if mod(quasar_ind,1000)==0
        fprintf('QSO %d of %d\n', num_quasars, num_quasars);
    end
end

            




