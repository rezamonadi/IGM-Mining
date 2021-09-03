% Build catalogs usable for spectra from dr16

Cooksey_C4_detected = fitsread(...
'data/C4_catalogs/Cooksey_C4_cat/distfiles/sdss_civ_cookseyetal13_update1.fit',...
'binarytable');

% Cooksey_C4_detected = fitsread(...
% '/home/reza/Desktop/sdss_civrate_hyb_all_update1.fit',...
% 'binarytable');
c4_QSO_ID                   = Cooksey_C4_detected{1};
% c4_detcted_mjd_dr7           = Cooksey_C4_detected{2};
% c4_detcted_plate_dr7         = Cooksey_C4_detected{3};
% c4_detcted_fiber_dr7         = Cooksey_C4_detected{4};
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
% all_snrs                = cooksey_catalog{172};
% all_bal_flag        = cooksey_catalog{131};
num_quasars             = numel(all_zqso);

all_z_c4 = zeros(num_quasars,1);
all_z_c4 = all_z_c4 -1;
all_NCIV = zeros(num_quasars,1);
all_c4_NCIV =zeros(num_quasars,1)-1;
all_QSO_ID=cell(num_quasars,1);
all_RATING=zeros(num_quasars,1)-1;
for i=1:num_quasars
    all_QSO_ID{i}=sprintf('%05i-%04i-%03i', (all_mjd_dr7(i)), ...
    (all_plate_dr7(i)), (all_fiber_dr7(i)));
end
%  adding a cloumn for c4 col density if there is a c4 for a sight line



% save catalog 
release = 'dr7';
variables_to_save = {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7', ...
 'all_QSO_ID', 'all_RA', 'all_DEC', 'all_zqso',...
  'EW1', 'EW2'};
save(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_save{:}, '-v7.3');

