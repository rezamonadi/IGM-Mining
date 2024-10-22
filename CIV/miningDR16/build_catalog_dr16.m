% Build catalogs usable for spectra from dr16

dr16q = ...
fitsread('data/DR16Q_v4.fits', 'binarytable');
all_plate_dr16             = dr16q{4};
all_mjd_dr16               = dr16q{5};
all_fiber_dr16             = dr16q{6};
all_RA_dr16                = dr16q{2};
all_DEC_dr16               = dr16q{3};
all_zqso_dr16              = dr16q{27};
all_BAL_PROB               = dr16q{57};



num_quasars_dr16 = numel(all_zqso_dr16);%num_quasars_dr16           = 1000;%numel(all_zqso_dr16);
all_QSO_ID_dr16=cell(num_quasars_dr16,1);


for i=1:num_quasars_dr16
    all_QSO_ID_dr16{i}=sprintf('%05i-%04i-%03i', (all_mjd_dr16(i)), ...
    (all_plate_dr16(i)), (all_fiber_dr16(i)));
end



% save catalog 
variables_to_save = {'all_plate_dr16', 'all_mjd_dr16', 'all_fiber_dr16', ...
 'all_QSO_ID_dr16', 'all_RA_dr16', 'all_DEC_dr16', 'all_zqso_dr16' };
save(sprintf('%s/catalog', processed_directory(releaseTest)), ...
    variables_to_save{:}, '-v7.3');
