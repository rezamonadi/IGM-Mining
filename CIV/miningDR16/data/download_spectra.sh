
# Download spectra of dr7 from the list provided by Cooksey et. al 
# http://www.guavanator.uhh.hawaii.edu/~kcooksey/SDSS/CIV/data/dr7qso_CIV_noBAL.list

pushd dr16

wget -x -nH -i ../dr16Spec.list 
# rsync -avzL --files-from=../dr16Spec.list rsync://dtn.sdss.org/dr16  . 2> /dev/null





