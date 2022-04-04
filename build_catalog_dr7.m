% Build catalogs usable for spectra from dr16

Cooksey_C4_detected = fitsread(...
'data/C4_catalogs/Cooksey_C4_cat/distfiles/sdss_civ_cookseyetal13_update1.fit',...
'binarytable');

% Cooksey_C4_detected = fitsread(...
% 'data/C4_catalogs/Cooksey_C4_cat/distfiles/sdss_civrate_hyb_all_update1.fit',...
% 'binarytable');


c4_QSO_ID                   = Cooksey_C4_detected{1};
% c4_detcted_mjd_dr7           = Cooksey_C4_detected{2};
% c4_detcted_plate_dr7         = Cooksey_C4_detected{3};
% c4_detcted_fiber_dr7         = Cooksey_C4_detected{4};
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
% ind = NCOLMFLG(:,1)==1 & NCOLMFLG(:,2)==1 & NCIV_ORG(:,1)>0 & NCIV_ORG(:,2)>0;

EW1                          = EW(:,1);
EW2                          = EW(:,2);
[nSys,dd]=size(c4_QSO_ID);
NCIV=zeros(nSys,1);
Z_c4=zeros(nSys,1);
for i=1:nSys
    NCIV(i) = NCIV_ORG(i,1)/SigmaNCIV_ORG(i,1)^2 + NCIV_ORG(i,2)/SigmaNCIV_ORG(i,2)^2;
    NCIV(i)=NCIV(i)/(1/SigmaNCIV_ORG(i,1)^2+1/SigmaNCIV_ORG(i,2)^2);
    Z_c4(i) = (Z_abs_ORG(i,1)+Z_abs_ORG(i,2))/2;
end




% f = fopen('data/C4_catalogs/Cooksey_C4_cat/processed/c4_catalog','w');
% for i=1:nSys
%       fprintf(f,'%05i-%04i-%03i  %f %f\n', c4_detcted_mjd_dr7(i), ...
%       c4_detcted_plate_dr7(i), c4_detcted_fiber_dr7(i), Z_c4(i), NCIV(i));
% end

save('data/C4_catalogs/Cooksey_C4_cat/processed/CIV-cat.mat','c4_QSO_ID','Z_c4','NCIV');

% There are some NAN valued c4_NCIV
% extract basic QSO information from Cookse_all_QSO catalog 
cooksey_catalog = ...
fitsread('data/dr7/distfiles/dr7qso_CIV_noBAL.fit', 'binarytable');
all_plate_dr7             = cooksey_catalog{48};
all_mjd_dr7             = cooksey_catalog{47};
all_fiber_dr7             = cooksey_catalog{49};
all_RA                = cooksey_catalog{2};
all_DEC               = cooksey_catalog{3};
all_zqso                = cooksey_catalog{4};
num_quasars             = numel(all_zqso);
all_QSO_ID=cell(num_quasars,1);

all_z_civ = zeros(num_quasars, 17)-1;
all_N_civ = zeros(num_quasars, 17)-1;
all_RATING = zeros(num_quasars, 17)-1;
for i=1:num_quasars
    all_QSO_ID{i}=sprintf('%05i-%04i-%03i', (all_mjd_dr7(i)), ...
    (all_plate_dr7(i)), (all_fiber_dr7(i)));
    ThisSystems = ismember(c4_QSO_ID, all_QSO_ID{i});
    thisZ_c4s = Z_c4(ThisSystems);
    thisN_c4s = NCIV(ThisSystems);
    this_RATING = RATING(ThisSystems);
    nSys = nnz(ThisSystems);
    
    for j=1:nSys
        all_z_civ(i,j) = thisZ_c4s(j);
        all_N_civ(i,j) = thisN_c4s(j);
        all_RATING(i,j) = this_RATING(j);
    end
end


% % dla catalog 
% dla_catalog = ...
% fitsread('data/dr7/distfiles/match-dla-civ.fits', 'binarytable');
% dla_plate  = dla_catalog{48};
% dla_mjd    = dla_catalog{47};
% dla_fiber  = dla_catalog{49};
% num_dlas   = length(dla_plate);
% log_posteriors_dla = dla_catalog{77};
% log_posteriors_no_dla = dla_catalog{78};
% dla_QSO_ID = cell(num_dlas,1);

% for i=1:num_dlas
%     dla_QSO_ID{i}=sprintf('%05i-%04i-%03i', (dla_mjd(i)), ...
%     (dla_plate(i)), (dla_fiber(i)));
% end


% save catalog 
release = 'dr7';
variables_to_save = {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7', ...
 'all_QSO_ID', 'all_RA', 'all_DEC', 'all_zqso', 'EW1', 'EW2',...
 'all_N_civ','all_z_civ', 'all_RATING', 'c4_QSO_ID'};
save(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_save{:}, '-v7.3');
