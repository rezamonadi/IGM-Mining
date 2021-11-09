% Build catalogs usable for spectra from dr16

dr12q = ...
fitsread('data/dr12/distfiles/dr12q_good.fits', 'binarytable');
all_plate_dr12             = dr12q{5};
all_mjd_dr12          = dr12q{6};
all_fiber_dr12             = dr12q{7};
all_RA_dr12                = dr12q{2};
all_DEC_dr12               = dr12q{3};
all_zqso_dr12                = dr12q{8};
num_quasars_dr12             = numel(all_zqso_dr12);
all_QSO_ID_dr12=cell(num_quasars_dr12,1);

for i=1:num_quasars_dr12
    all_QSO_ID{i}=sprintf('%05i-%04i-%03i', (all_mjd_dr12(i)), ...
    (all_plate_dr12(i)), (all_fiber_dr12(i)));
end



% save catalog 
variables_to_save = {'all_plate_dr12', 'all_mjd_dr12', 'all_fiber_dr12', ...
 'all_QSO_ID_dr12', 'all_RA_dr12', 'all_DEC_dr12', 'all_zqso_dr12' };
save(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_save{:}, '-v7.3');
