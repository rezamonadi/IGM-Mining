% Build catalogs usable for spectra from dr16

Cooksey_C4_detected = fitsread(...
'data/C4_catalogs/Cooksey_C4_cat/distfiles/sdss_mgii_seyffertetal13.fit',...
'binarytable');
c4_QSO_ID                    = Cooksey_C4_detected{1};
z_qso_system                 = Cooksey_C4_detected{10};
Z_abs_ORG                    = Cooksey_C4_detected{17};
EW                           = Cooksey_C4_detected{22};
SigmaEW                      = Cooksey_C4_detected{23};
flagEW                       = Cooksey_C4_detected{24};
NCIV_ORG                     = Cooksey_C4_detected{27};
SigmaNCIV_ORG                = Cooksey_C4_detected{28};
NCOLMFLG                     = Cooksey_C4_detected{29};
dummy                        = Cooksey_C4_detected{30};

RATING                       = dummy(:,1);

% % filtering out those column densities with not good measurements
EW1                          = EW(:,1);
errEW1                       = SigmaEW(:,1);
errEW2                     = SigmaEW(:,2);
EW2                          = EW(:,2);
[nSys,dd]=size(c4_QSO_ID);
NCIV=zeros(nSys,1);
Z_c4=zeros(nSys,3);
for i=1:nSys
    NCIV(i) = NCIV_ORG(i,1)/SigmaNCIV_ORG(i,1)^2 + NCIV_ORG(i,2)/SigmaNCIV_ORG(i,2)^2;
    NCIV(i)=NCIV(i)/(1/SigmaNCIV_ORG(i,1)^2+1/SigmaNCIV_ORG(i,2)^2);
    Z_c4(i,1) = min([Z_abs_ORG(i,1), Z_abs_ORG(i,2)]);
    Z_c4(i,2) = max([Z_abs_ORG(i,1), Z_abs_ORG(i,2)]);
    Z_c4(i,3) = mean([Z_abs_ORG(i,1), Z_abs_ORG(i,2)]);
end

save('data/C4_catalogs/Cooksey_C4_cat/processed/CIV-cat.mat','c4_QSO_ID','Z_c4','NCIV');

% There are some NAN valued c4_NCIV
% extract basic QSO information from Cookse_all_QSO catalog 
cooksey_catalog = ...
fitsread('data/dr7/distfiles/dr7_QSO_MgII.fits', 'binarytable');
all_plate_dr7             = cooksey_catalog{48};
all_mjd_dr7             = cooksey_catalog{47};
all_fiber_dr7             = cooksey_catalog{49};
all_RA                = cooksey_catalog{2};
all_DEC               = cooksey_catalog{3};
all_zqso                = cooksey_catalog{4};
num_quasars             = numel(all_zqso);
all_QSO_ID=cell(num_quasars,1);

all_z_civ1 = zeros(num_quasars, 17)-1;
all_z_civ2 = zeros(num_quasars, 17)-1;
all_z_civ3 = zeros(num_quasars, 17)-1;
all_N_civ = zeros(num_quasars, 17)-1;
all_RATING = zeros(num_quasars, 17)-1;
all_EW1 = zeros(num_quasars, 17)-1;
all_EW2 = zeros(num_quasars, 17)-1;
all_errEW1 = zeros(num_quasars,17)-1;
all_errEW2 = zeros(num_quasars,17)-1;
for i=1:num_quasars
    all_QSO_ID{i}=sprintf('%05i-%04i-%03i', (all_mjd_dr7(i)), ...
    (all_plate_dr7(i)), (all_fiber_dr7(i)));
    ThisSystems = ismember(c4_QSO_ID, all_QSO_ID{i});
    thisZ_c4s_1 = Z_c4(ThisSystems,1);
    thisZ_c4s_2 = Z_c4(ThisSystems,2);
    thisZ_c4s_3 = Z_c4(ThisSystems,3);
    thisN_c4s = NCIV(ThisSystems);
    this_RATING = RATING(ThisSystems);
    this_EW1 = EW1(ThisSystems);
    this_EW2 = EW1(ThisSystems);
    this_errEW1 = errEW1(ThisSystems);
    this_errEW2 = errEW2(ThisSystems);
    nSys = nnz(ThisSystems);
    
    for j=1:nSys
        all_z_civ1(i,j) = thisZ_c4s_1(j);
        all_z_civ2(i,j) = thisZ_c4s_2(j);
        all_z_civ3(i,j) = thisZ_c4s_3(j);
        all_N_civ(i,j) = thisN_c4s(j);
        all_RATING(i,j) = this_RATING(j);
        all_EW1(i, j) = this_EW1(j);
        all_EW2(i, j) = this_EW2(j);
        all_errEW1(i,j)= this_errEW1(j);
        all_errEW2(i,j)= this_errEW2(j);

    end
end



% save catalog 
release = 'dr7';
variables_to_save = {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7', ...
 'all_QSO_ID', 'all_RA', 'all_DEC', 'all_zqso', 'all_EW1', 'all_EW2', 'all_errEW1', 'all_errEW2', ...
 'all_N_civ','all_z_civ1', 'all_z_civ2', 'all_z_civ3' 'all_RATING', 'c4_QSO_ID'};
save(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_save{:}, '-v7.3');
