# clean.bash, filters, runs deflare and generate exposure map

# Filter    
dmcopy "*_repro_evt2.fits[energy=300:9000,ccd_id=$2]" $1_evt2_c$2.fits clobber=yes

# Extract lightcurve and filter with deflare
dmextract "$1_evt2_c$2.fits[bin time=::200]" outfile=lc.fits opt=ltc1
deflare "lc.fits[count_rate<2]" outfile=lc_2.gti method=clean
dmcopy $1_evt2_c$2.fits evt2_filt.fits
dmcopy "$1_evt2_c$2.fits[@lc_2.gti]" evt2_filt.fits clobber=yes
echo "==========================================================
Ignore any errors here, indicates an acceptable lightcurve
=========================================================="

# Create exposure map and image
fluximage evt2_filt.fits binsize=1 bands=broad outroot=image psfecf=0.393 clobber=yes