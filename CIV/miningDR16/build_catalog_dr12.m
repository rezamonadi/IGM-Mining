% Build catalogs usable for spectra from dr16

dr12q = ...
fitsread('data/dr16/distfiles/dr16_QSO_noBAL.fit', 'binarytable');
all_plate_dr12             = dr12q{4};
all_mjd_dr12               = dr12q{5};
all_fiber_dr12             = dr12q{6};
all_RA_dr12                = dr12q{2};
all_DEC_dr12               = dr12q{3};
all_zqso_dr12              = dr12q{27};
num_quasars_dr12 = min(2000, numel(all_zqso_dr12));%num_quasars_dr12           = 1000;%numel(all_zqso_dr12);
all_QSO_ID_dr12=cell(num_quasars_dr12,1);

subset_all_zqso_dr12 = all_zqso_dr12(1:num_quasars_dr12);
all_zqso_dr12 = subset_all_zqso_dr12; 
% Extract the first 1000 entries (or fewer if less than 1000 entries exist)

for i=1:num_quasars_dr12
    all_QSO_ID_dr12{i}=sprintf('%05i-%04i-%03i', (all_mjd_dr12(i)), ...
    (all_plate_dr12(i)), (all_fiber_dr12(i)));
end



% save catalog 
variables_to_save = {'all_plate_dr12', 'all_mjd_dr12', 'all_fiber_dr12', ...
 'all_QSO_ID_dr12', 'all_RA_dr12', 'all_DEC_dr12', 'all_zqso_dr12' };
save(sprintf('%s/catalog', processed_directory(releaseTest)), ...
    variables_to_save{:}, '-v7.3');
