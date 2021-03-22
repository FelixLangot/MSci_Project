# This program extracts a spectrum from an observation

echo "Do you want to download and reprocess data? ('y' for yes, 'n' for no)"
read ans 

if [ "$ans" = "y" ]; then
    # Downloads argument obsid
    download_chandra_obsid $1

    #moves into directory and reprocesses
    cd $1/
    chandra_repro indir=. outdir=./repro/
else
cd $1/

fi

# Navigates into directory and reprocesses
cd repro/
punlearn dmcopy

## Open ds9 and choose ccd
#ds9 *_repro_evt2* &
#echo "Which ccd to use ?"   
#read ccd
    
dmcopy "*_repro_evt2.fits[energy=500:2000]" $1_evt2.fits clobber=yes

# Open ds9 and create src and bkg_2 regions
ds9 $1_evt2.fits &

read -p "Construct circular src.reg and bkg_2.reg and press any key to continue" -n1 -s
echo

# Option to remove annuli
echo "Do you wish to remove a point source in annuli ? ('y' or 'n')"
read rem
if [ "$rem" = "y" ]; then 
    ds9 $1_evt2.fits &
    read -p "Construct a region to remove and press any key to continue. (name = contam.reg)" -n1 -s
    echo
    dmcopy "$1_evt2.fits[exclude sky=region(contam.reg)]" $1_evt2.fits clobber=yes
elif [ "$rem" = "n" ]; then
    dmcopy "$1_evt2.fits" $1_evt2.fits clobber=yes
fi
echo "Proceding to run specextract ... "
sleep 1.5
# Run specextract
specextract infile="$1_evt2.fits[sky=region(src.reg)]" outroot=spec bkgfile="$1_evt2.fits[sky=region(bkg_2.reg)]" clobber=yes

echo "specextract done"

# Run Sherpa script
echo "Now run this Sherpa script with determined redshift and nH values:"
echo "SHERPA SCRIPT - TEMPERATURE =============================================
import matplotlib.pyplot as plt
load_pha("\""spec.pi"\"")
notice(0.5,2)
set_ylog()
set_source(xsphabs.abs1 * xsapec.apec)
thaw(apec.Abundanc)
abs1.nH = CHANGEME
apec.redshift = CHANGEME
freeze(abs1.nH)
freeze(apec.redshift)
guess(apec)
set_stat("\""wstat"\"")
fit()
covar()
x = get_covar_results()
plot_fit()
plt.savefig("\""../../$1T_eprofile.eps"\"", dpi=300)
subtract()
plot_fit()
plt.savefig("\""../../$1T_eprofilesub.eps"\"", dpi=300)
============================================="
echo "SHERPA SCRIPT - COOLING FUNCTION =============================================       
apec.redshift = 0   # for cooling constant
apec.norm = 1
F = calc_photon_flux(0.5 * (1+z[CHANGEME]), 2 * (1+z[CHANGEME]))   # * (1+z) to find in cluster frame
eF = F*x.parmaxes[2]/x.parvals[2]
print('flux F = ', F, '±', eF, 'photons/cm^2/s')
Λ = F * 10**(-14)
eΛ = eF * 10**(-14)
print('cooling function Λ = ', Λ, '±', eΛ)
============================================="
sherpa 