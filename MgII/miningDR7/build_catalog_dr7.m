% Build catalogs usable for spectra from dr16

Seyfert_MgII_detected = fitsread(...
'data/MgII_catalogs/Seyfert_MgII_cat/distfiles/sdss_mgii_seyffertetal13.fit',...
'binarytable');
MgII_QSO_ID                  = Seyfert_MgII_detected{1};
z_qso_system                 = Seyfert_MgII_detected{10};
Z_abs_ORG                    = Seyfert_MgII_detected{17};
EW                           = Seyfert_MgII_detected{22};
SigmaEW                      = Seyfert_MgII_detected{23};
flagEW                       = Seyfert_MgII_detected{24};
NMgII_ORG                    = Seyfert_MgII_detected{27};
SigmaNMgII_ORG               = Seyfert_MgII_detected{28};
NCOLMFLG                     = Seyfert_MgII_detected{29};
dummy                        = Seyfert_MgII_detected{30};

RATING                       = dummy(:,1);

% % filtering out those column densities with not good measurements
EW1                          = EW(:,1);
errEW1                       = SigmaEW(:,1);
errEW2                     = SigmaEW(:,2);
EW2                          = EW(:,2);
[nSys,dd]=size(MgII_QSO_ID);
NMgII=zeros(nSys,1);
Z_MgII=zeros(nSys,3);
for i=1:nSys
    NMgII(i) = NMgII_ORG(i,1)/SigmaNMgII_ORG(i,1)^2 + NMgII_ORG(i,2)/SigmaNMgII_ORG(i,2)^2;
    NMgII(i)= NMgII(i)/(1/SigmaNMgII_ORG(i,1)^2+1/SigmaNMgII_ORG(i,2)^2);
    Z_MgII(i,1) = min([Z_abs_ORG(i,1), Z_abs_ORG(i,2)]);
    Z_MgII(i,2) = max([Z_abs_ORG(i,1), Z_abs_ORG(i,2)]);
    Z_MgII(i,3) = mean([Z_abs_ORG(i,1), Z_abs_ORG(i,2)]);
end

save('data/MgII_catalogs/Seyfert_MgII_cat/processed/MgII-cat.mat','MgII_QSO_ID','Z_MgII','NMgII');

% There are some NAN valued c4_NCIV
% extract basic QSO information from Cookse_all_QSO catalog 
QSO_catalog = ...
fitsread('data/dr7/distfiles/dr7_QSO_MgII.fits', 'binarytable');
all_plate_dr7             = QSO_catalog{48};
all_mjd_dr7               = QSO_catalog{47};
all_fiber_dr7             = QSO_catalog{49};
all_RA                    = QSO_catalog{2};
all_DEC                   = QSO_catalog{3};
all_zqso                  = QSO_catalog{4};
num_quasars               = numel(all_zqso);
all_QSO_ID                = cell(num_quasars,1);
all_z_MgII1               = zeros(num_quasars, 17)-1;
all_z_MgII2               = zeros(num_quasars, 17)-1;
all_z_MgII3               = zeros(num_quasars, 17)-1;        
all_MgII2                 = zeros(num_quasars, 17)-1;
all_MgII3                 = zeros(num_quasars, 17)-1;
all_N_MgII                = zeros(num_quasars, 17)-1;
all_RATING                = zeros(num_quasars, 17)-1;
all_EW1                   = zeros(num_quasars, 17)-1;
all_EW2                   = zeros(num_quasars, 17)-1;
all_errEW1                = zeros(num_quasars,17)-1;
all_errEW2                = zeros(num_quasars,17)-1;
for i=1:num_quasars
    all_QSO_ID{i}=sprintf('%05i-%04i-%03i', (all_mjd_dr7(i)), ...
    (all_plate_dr7(i)), (all_fiber_dr7(i)));
    ThisSystems = ismember(MgII_QSO_ID, all_QSO_ID{i});
    thisZ_MgIIs_1 = Z_MgII(ThisSystems,1);
    thisZ_MgIIs_2 = Z_MgII(ThisSystems,2);
    thisZ_MgIIs_3 = Z_MgII(ThisSystems,3);
    thisN_MgIIs = NMgII(ThisSystems);
    this_RATING = RATING(ThisSystems);
    this_EW1 = EW1(ThisSystems);
    this_EW2 = EW1(ThisSystems);
    this_errEW1 = errEW1(ThisSystems);
    this_errEW2 = errEW2(ThisSystems);
    nSys = nnz(ThisSystems);
    
    for j=1:nSys
        all_z_MgII1(i,j) = thisZ_MgIIs_1(j);
        all_z_MgII2(i,j) = thisZ_MgIIs_2(j);
        all_z_MgII3(i,j) = thisZ_MgIIs_3(j);
        all_N_MgII(i,j) = thisN_MgIIs(j);
        all_RATING(i,j) = this_RATING(j);
        all_EW1(i, j) = this_EW1(j);
        all_EW2(i, j) = this_EW2(j);
        all_errEW1(i,j)= this_errEW1(j);
        all_errEW2(i,j)= this_errEW2(j);

    end
end



% save catalog 
release = 'dr7';
variables_to_save = {'all_plate_dr7',...
                     'all_mjd_dr7',...
                     'all_fiber_dr7', ...
                     'all_QSO_ID', ...
                     'all_RA', ...
                     'all_DEC', ...
                     'all_zqso', ...
                     'all_EW1', ...
                     'all_EW2', ...
                     'all_errEW1', ...
                     'all_errEW2', ...
                     'all_N_MgII', ...
                     'all_z_MgII1', ...
                     'all_z_MgII2', ...
                     'all_z_MgII3', ...
                     'all_RATING', ...
                     'MgII_QSO_ID'};
                     
save(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_save{:}, '-v7.3');
