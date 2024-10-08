readin = 0;
if readin == 1
data = fitsread('dr16_QSO_noBAL.fit', 'binarytable');
end
plate = data{:,4};
MJD = data{:,5};
FID = data{:,6};

fileIDI = fopen('MJDPLATEFID.txt','w');
astrotable = table(MJD, plate, FID);
writetable(astrotable, 'MJDPLATEFID.txt');

fileIDII = fopen('dr16_CIV_noBAL.list', 'w');
for i = 1:numel(plate)
    file_name = sprintf('https://dr16.sdss.org/sas/dr16/sdss/spectro/redux/v5_13_0/spectra/lite/%d/spec-%d-%d-%.4d.fits\n', plate(i),plate(i),MJD(i), FID(i));
    fprintf(fileIDII, file_name);
end