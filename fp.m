% False plositve plotter 

clear
set_parameters;
build_catalog;

% training_set_name ='-WNciv-to-1600';
filename = sprintf('%s/processed_qsos_R%s', ...
    processed_directory(release), ...
    training_set_name);
load(filename, 'test_ind', 'p_L1','p_no_c4', 'map_z_c4L2', 'training_set_name',...
                'map_N_c4L2', 'map_sigma_c4L2');
% training_set_name
cooksey_catalog = ...
fitsread('data/dr7/distfiles/dr7qso_CIV_noBAL.fit', 'binarytable');
all_zqso          = cooksey_catalog{4};
ind_has_c4 =ismember(all_QSO_ID, c4_QSO_ID) & (filter_flags==0) & (all_RATING==3);

a = 1e-5;
p_c4 = 1- p_no_c4 -a*p_L1;
tr=0.9;
min_N=linspace(min(map_N_c4L2), max(map_N_c4L2)-0.1, 100);
FP = zeros(100,1);
for i=1:100

    FP(i) = nnz(~ind_has_c4(test_ind) & p_c4>tr & map_N_c4L2>min_N(i) );
end
figure()
plot(min_N, FP)
set(get(gca, 'XLabel'), 'String', 'min(N)');
set(get(gca, 'YLabel'), 'String', 'FP');
saveas(gca, 'minN-FP.png')


Dz=linspace(0, 0.2, 100);
FP = zeros(100,1);
for i=1:100

    FP(i) = nnz(~ind_has_c4(test_ind) & p_c4>tr & (all_zqso(test_ind)> (map_z_c4L2+Dz(i))));
end
figure()
plot(Dz, FP)
set(get(gca, 'XLabel'), 'String', 'Dz');
set(get(gca, 'YLabel'), 'String', 'FP');
saveas(gca, 'Dz-FP.png')


min_sigma=linspace(min(map_sigma_c4L2), max(map_sigma_c4L2)-5, 100);
FP = zeros(100,1);
for i=1:100

    FP(i) = nnz(~ind_has_c4(test_ind) & p_c4>tr & map_sigma_c4L2>min_sigma(i) );
end
figure()
plot(min_sigma, FP)
set(get(gca, 'XLabel'), 'String', 'min(\sigma)');
set(get(gca, 'YLabel'), 'String', 'FP');
saveas(gca, 'minSigma-FP.png')
