# This program extracts a radial profile from an observation

# Download data?
bash ~/Data/dl.bash $1

cd $1/repro/

# Open ds9 and choose ccd
ds9 *_repro_evt2* &
echo "Which ccd to use?" 
read ccd

# Filter, run deflare and generate exposure map
echo "Do you want to clean data ('y' or 'n')"
read cl
if [ "$cl" = "y" ]; then 
    bash ~/Data/clean.bash $1 $ccd
fi

# Open ds9 and create annuli and bkg regions
ds9 evt2_filt.fits &

read -p "Construct annuli.reg and bkg.reg and press any key to continue" -n1 -s
echo

# Option to remove annuli
echo "Do you wish to remove a point source in annuli ? ('y' or 'n')"
read rem
if [ "$rem" = "y" ]; then 
    ds9 evt2_filt.fits &
    read -p "Construct a region to remove and press any key to continue. (name = contam.reg)" -n1 -s
    echo
    dmcopy "evt2_filt.fits[exclude sky=region(contam.reg)]" $1_evt2_ready_c$ccd.fits clobber=yes
else
    dmcopy "evt2_filt.fits" $1_evt2_ready_c$ccd.fits clobber=yes
fi

# Run dmextract
dmextract infile="$1_evt2_ready_c$ccd.fits[bin sky=@annuli.reg]" outfile=$1_rprofile.fits bkg="$1_evt2_ready_c$ccd.fits[bin sky=@bkg.reg]" exp=image_broad_thresh.expmap bkgexp=image_broad_thresh.expmap opt=generic clobber=yes
dmextract infile="image_broad_thresh.expmap[bin sky=@annuli.reg]" outfile=exp_rprofile.fits bkg="image_broad_thresh.expmap[bin sky=@bkg.reg]" opt=generic clobber=yes exp=image_broad_thresh.expmap bkgexp=image_broad_thresh.expmap clobber=yes
echo "==================
dmextract finished
=================="

# Run Sherpa script
bash ~/Data/sherps.bash $1

