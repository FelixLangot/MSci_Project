################################## This program extracts a spectrum from an observation ##################################

# Navigates into directory
cd $1/repro/

echo "Do you want to draw regions? ('y' or 'n')" 
read draw
if [ "$draw" = "y" ]; then
    punlearn dmcopy
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
elif [ "$draw" = "n" ]; then
    echo "Using already drawn regions"
    find contam.reg
fi

echo 'Run specextract ? (y or n) '
read spec
if [ "$spec" = "y" ]; then
    echo "Proceding to run specextract ... "
    sleep 1
    # Run specextract
    specextract infile="*_repro_evt2.fits[sky=region(src.reg)]" outroot=spec bkgfile="*_repro_evt2.fits[sky=region(bkg_2.reg)]" clobber=yes
    echo "specextract done"
elif [ "$spec" = "n" ]; then
    echo "Proceding to graphing and analysis"
    sleep 1
fi

# Run Sherpa script
echo "Analysisng cluster in obsID '$1'"
echo "nH ?"
read nHc
echo "z ?"
read zred
echo "
import time
print('====================================== SHERPA SCRIPT - TEMPERATURE AND COOLING FUNCTION ======================================')
import matplotlib.pyplot as plt
print('---------------- TEMPERATURE ----------------')
load_pha("\""spec.pi"\"")
notice(0.5,7)
set_ylog()
set_source(xsphabs.abs1 * xsapec.apec)
thaw(apec.Abundanc)
abs1.nH = $nHc
apec.redshift = $zred
z = apec.redshift
freeze(abs1.nH)
freeze(apec.redshift)
guess(apec)
set_stat("\""wstat"\"")
fit()
covar()
x = get_covar_results()
b = 30
group_counts(1, b)
plot_fit()
plt.savefig("\""../$1T_eprofile_"\""+str(b).zfill(2)+"\""cntspbn.eps"\"", dpi=300)
subtract()
plot_fit()
plt.title('$1')
plt.savefig("\""../$1T_eprofilesub_"\""+str(b).zfill(2)+"\""cntspbn.eps"\"", dpi=300)  
print('Plotted with ', b, ' counts per bin.') 

plotinput = '0'
while plotinput != 'q':
    print('If you want to change the number of counts per bin, type a new number (last = ', b,'). If you do not want to change, type q. ')
    plotinput = input()
    if plotinput != 'q':
        b = int(plotinput)
        group_counts(1, b)
        unsubtract()
        plot_fit()
        plt.title('$1')
        plt.savefig("\""../$1T_eprofile_"\""+str(b).zfill(2)+"\""cntspbn.eps"\"", dpi=300)
        subtract()
        plot_fit()
        plt.title('$1')
        plt.savefig("\""../$1T_eprofilesub_"\""+str(b).zfill(2)+"\""cntspbn.eps"\"", dpi=300)  
        print('Plotted with ', b, ' counts per bin.') 
    elif plotinput == 'q':
        print('Proceeding to cooling function')
        time.sleep(1)

print('---------------- COOLING FUNCTION ----------------')
zsaved = z.val
apec.redshift = 0   # for cooling constant
apec.norm = 1
F = calc_photon_flux(0.5 * (1+zsaved), 2 * (1+zsaved))   # * (1+z) to find in cluster frame
eF = F*x.parmaxes[2]/x.parvals[2]
Λ = F * 10**(-14)
eΛ = eF * 10**(-14)
print('Λ = ', Λ, '±', eΛ)
print('ObsId = $1')
print('==============================================================================================')
exit
" > sherpaTscript.py
sherpa -n sherpaTscript.py 