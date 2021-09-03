clear
clc 
set_parameters
build_catalog
training_set_name = 'Rvoigt-10000Smp0-tr-80-20-vCutKathy-3000';
filename = sprintf('%s/processed_qsos_R%s.mat', ...
    processed_directory(release), ...
    training_set_name);
load(filename);
mkdir(sprintf('ROC-R-%s', training_set_name));
ind_has_c4 =ismember(all_QSO_ID, c4_QSO_ID);

for a =[0, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]
    fig= figure();

    p_c4New = 1- p_no_c4 - a*p_L1;
    y_score = p_c4New;
    y_true = ind_has_c4(test_ind);
    [X,Y,T,AUC] =perfcurve(y_true, y_score, 'true');
    plot(X,Y)
    legend(sprintf('a=%.2e, AUC=%.5f',a, AUC))
    set(get(gca, 'YLabel'), 'String', 'TPR');
    set(get(gca, 'XLabel'), 'String', 'FPR');
    exportgraphics(fig, sprintf('ROC-R-%s/a-%.2e.pdf', training_set_name, a)...
                  ,'ContentType','vector')
end
