% Build catalogs usable for spectra from dr16

dr16q = ...
fitsread('data/dr16/distfiles/dr16_QSO_noBAL.fit', 'binarytable');
all_plate_dr16             = dr16q{4};
all_mjd_dr16               = dr16q{5};
all_fiber_dr16             = dr16q{6};
all_RA_dr16                = dr16q{2};
all_DEC_dr16               = dr16q{3};
all_zqso_dr16              = dr16q{27};

% Extract the first 1000 entries (or fewer if less than 1000 entries exist)
ind = [];
for i = 1:numel(all_plate_dr16)
    if all_plate_dr16(i)>=3586
        ind(end+1) = i;
    end
end
all_plate_dr16 = all_plate_dr16(ind);
all_mjd_dr16 = all_mjd_dr16(ind);
all_fiber_dr16 = all_fiber_dr16(ind);
all_RA_dr16 = all_RA_dr16(ind);
all_DEC_dr16 = all_DEC_dr16(ind);
all_zqso_dr16 = all_zqso_dr16(ind);

num_quasars_dr16 = min(2000, numel(all_zqso_dr16));%num_quasars_dr16           = 1000;%numel(all_zqso_dr16);
all_QSO_ID_dr16=cell(num_quasars_dr16,1);

subset_all_zqso_dr16 = all_zqso_dr16(1:num_quasars_dr16);
all_zqso_dr16 = subset_all_zqso_dr16; 

for i=1:num_quasars_dr16
    all_QSO_ID_dr16{i}=sprintf('%05i-%04i-%03i', (all_mjd_dr16(i)), ...
    (all_plate_dr16(i)), (all_fiber_dr16(i)));
end



% save catalog 
variables_to_save = {'all_plate_dr16', 'all_mjd_dr16', 'all_fiber_dr16', ...
 'all_QSO_ID_dr16', 'all_RA_dr16', 'all_DEC_dr16', 'all_zqso_dr16' };
save(sprintf('%s/catalog', processed_directory(releaseTest)), ...
    variables_to_save{:}, '-v7.3');
