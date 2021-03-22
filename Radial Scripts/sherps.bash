#sherps.bash, this prints a sherpa script

echo "Method 1 (Sherpa emap), method 2 (manual emap) or method 3 (no emap)? ('1', '2' or '3')"
read method

if [ "$method" = "1" ]; then
    echo "load_table_model(\"emap\", \"exp_rprofile.fits[cols SUR_BRI]\")
load_data(1,"\""$1_rprofile.fits"\"", 3, ["\""RMID"\"","\""SUR_BRI"\"","\""SUR_BRI_ERR"\""])
set_source(beta1d.src * emap)
src.r0=10
src.beta=0.6
src.ampl=1
freeze(emap.ampl)
set_xlog()
set_ylog()
fit()
plot_fit()
" > sherpa_script.py

elif [ "$method" = "2" ]; then
    # Isolates exposure sur_bri
    dmcopy "exp_rprofile.fits[cols exp_sb=sur_bri]" renamed.fits clobber=yes

    # Appends sur_bri column to regular rprofile
    dmpaste infile=$1_rprofile.fits pastefile=renamed.fits outfile=combined.fits clobber=yes

    # Divides sur_bri and sur_bri_err by exposure amount
    dmtcalc infile=combined.fits outfile=almost.fits expression="corr=SUR_BRI/exp_sb" clobber=yes
    dmtcalc infile=almost.fits outfile=exp_corrected.fits expression="corr_err=SUR_BRI_ERR/exp_sb" clobber=yes

    echo "load_data(1,"\""exp_corrected.fits"\"", 3, ["\""RMID"\"","\""CORR"\"","\""CORR_ERR"\""])
set_source(beta1d.src)
src.r0=10
src.beta=0.6
src.ampl=1
set_xlog()
set_ylog()
fit()
plot_fit()
" > sherpa_script.py

else
    echo "load_data(1,"\""$1_rprofile.fits"\"", 3, ["\""RMID"\"","\""SUR_FLUX"\"","\""SUR_FLUX_ERR"\""])
set_source(beta1d.src)
src.r0=5
src.beta=0.6
src.ampl=0.1
set_xlog()
set_ylog()
fit()
plot_fit()
" > sherpa_script.py

fi

sherpa sherpa_script.py
