clear
clc 
set_parameters
build_catalog
% training_set_name = ''
filename = sprintf('%s/processed_qsos_R%s.mat', ...
    processed_directory(release), ...
    training_set_name);
load(filename);
Dz= 0.1
mkdir(sprintf('ROC-R-dz-%d-%s',fix(100*Dz), training_set_name));

mask_dz =(all_zqso(test_ind)> (map_z_c4L2+Dz));
for a =[0, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]
    fig= figure();

    p_c4New = 1- p_no_c4 - a*p_L1;
    y_score = p_c4New(mask_dz);
    y_true = all_ind_c4(test_ind);
    y_true = y_true(all_zqso(test_ind)> (map_z_c4L2+Dz));
    [X,Y,T,AUC] =perfcurve(y_true, y_score, 'true');
    plot(X,Y)
    legend(sprintf('a=%.2e, AUC=%.5f',a, AUC))
    set(get(gca, 'YLabel'), 'String', 'TPR');
    set(get(gca, 'XLabel'), 'String', 'FPR');
    exportgraphics(fig, sprintf('ROC-R-dz-%d-%s/a-%.2e.pdf',fix(Dz*100), training_set_name, a)...
                  ,'ContentType','vector')
end
