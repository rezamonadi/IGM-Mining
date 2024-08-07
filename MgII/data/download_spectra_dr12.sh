pushd spectra
rsync --info=progress2 -h --no-motd --files-from=dr12SpecList rsync://data.sdss.org/dr12/boss/spectro/redux/ . 2> /dev/null
