ObsIDarray=(13405 13465 13403 13396 13463 13402 13486 13398 13470 13487 13466 13476 13474 13488 14448 13397 13483 13491 13472 15573 13490 13485 13478 13467 13471 13473 13479 13493 14380 13489)

echo "All cluster IDs:"
echo ${ObsIDarray[*]}

echo "--------------------------------------------------------------------------"

echo "Total number of clusters:"
echo "${#ObsIDarray[@]}"
TotN=${#ObsIDarray[@]}
echo $TotN

echo "ObsIDarray[0]:"
echo "${ObsIDarray[0]}"

mkdir 'Clusters'

cd Clusters

for (( i=0; i<=$TotN; i++ ))
do
	k=${ObsIDarray[i]}
	echo -n "$i "
	download_chandra_obsid $k
	# moves into directory and reprocesses
    cd $k
    chandra_repro indir=. outdir=./repro/
	cd ..
done


#for i in range TotN
#	download_chandra_obsid iname
#	#moves into directory and reprocesses
#    cd iname
#    chandra_repro indir=. outdir=./repro/
