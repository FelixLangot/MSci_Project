# dl.bash, downloads and repros data

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

# Navigates into directory
cd repro/
punlearn dmcopy
