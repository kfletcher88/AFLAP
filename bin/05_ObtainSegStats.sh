#!/bin/bash -l
#################################################
#       A shell script to obtain segregation statistics and exclude progeny which have low coverage
#################################################
Name=$(basename $0)
usage="${Name}; [-h] [-P] [-m] -- A script to plot marker distributions in progeny (AFLAP 5/6).
Options
        -h show this help message
        -P Pedigree file, required. See AFLAP README for more information.
        -m K-mer size. Optional. Default [31]
Temporary files will be output to AFLAP_tmp/05."
#Option block
while getopts ':LhP:m:d:D:' option; do
        case "$option" in
                h)  echo "$usage"
                         exit
                         ;;
                P)  Ped=$OPTARG
                         ;;
                m)  mer=$OPTARG
                         ;;
		L)  CovCut=1
		    echo "Low Coverage setting used"
			 ;;
		D)  SegDistU=$OPTARG
			 ;;
		d)  SegDistL=$OPTARG
			 ;;
                \?) printf "illegal option: -%s\n\n" "$OPTARG" >&2
                    echo "$usage"
                    exit 1
                         ;;
        esac
done

#SetDefault CovCut
if [[ -z $CovCut ]] ; then CovCut=2 ; echo -e "\nDefault coverage cut-off = 2." ; fi
if [[ -e AFLAP_tmp/01/LA.txt && -e AFLAP_tmp/02/Boundaries.txt ]]
then
echo -e "\nIntermediate files detected"
else
echo -e "\n Could not find output of previous scripts. Please rerun full pipeline."
exit 1
fi
mkdir -p AFLAP_tmp/05/FilteredCall
mkdir -p AFLAP_tmp/05/SegregationInformation
mkdir -p AFLAP_Results
#BASH math
kk=$(printf %.0f $(echo $mer+$mer-1 | bc -l))

#Find bin directory for Rscripts
DIR=$(dirname $0)

#Loop through LA to run iteration for each parent required
for g in `cat AFLAP_tmp/01/LA.txt`
do
#1 Check for Genotype table
#Again need a low vs high coverage flag
	P0=$(cat AFLAP_tmp/03/${g}_CrossedTo.txt | tr '\n' '_' | sed 's/_$//')
	Lo=$(awk -v var="$g" '$1 == var {print $2}' AFLAP_tmp/02/Boundaries.txt)
	Up=$(awk -v var="$g" '$1 == var {print $3}' AFLAP_tmp/02/Boundaries.txt)
	if [[ -e AFLAP_tmp/04/${g}_m${mer}_L${Lo}_U${Up}_$P0.Genotypes.MarkerID.tsv ]]
	then
#Progeny counter
		ProC=$(awk 'NR == 1 {print NF-2}' AFLAP_tmp/04/${g}_m${mer}_L${Lo}_U${Up}_$P0.Genotypes.MarkerID.tsv)
		echo "Genotype calls for $g detected. Summarizing"
		sed 's/_/ /' AFLAP_tmp/04/${g}_m${mer}_L${Lo}_U${Up}_$P0.Genotypes.MarkerID.tsv | awk -v var="$kk" '$3 == var {for (i=4; i<=NF;i++) j+=$i; print j; j=0 }' | sort -n | uniq -c | awk -v var=$ProC '{print $2/var, $1}' > AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerEqual${kk}.hist
		sed 's/_/ /' AFLAP_tmp/04/${g}_m${mer}_L${Lo}_U${Up}_$P0.Genotypes.MarkerID.tsv | awk -v var="$kk" '$3 > var {for (i=4; i<=NF;i++) j+=$i; print j; j=0 }' | sort -n | uniq -c | awk -v var=$ProC '{print $2/var, $1}' > AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerOver${kk}.hist
		sed 's/_/ /' AFLAP_tmp/04/${g}_m${mer}_L${Lo}_U${Up}_$P0.Genotypes.MarkerID.tsv | awk '{for (i=4; i<=NF;i++) j+=$i; print j; j=0 }' | sort -n | uniq -c | awk -v var=$ProC '{print $2/var, $1}' > AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_AllMarkers.hist
		Rscript $DIR/SegStats.R AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerEqual${kk}.hist AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerOver${kk}.hist  AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_AllMarkers.hist AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerSeg.png
		for h in `awk -v var="$g" '$4 == var || $5 == var {print $1}' $Ped | sort -u`
		do
			if [[ -e AFLAP_tmp/04/Call/${h}_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt ]]
			then
			Mcou=$(awk '{sum += $1}END{print sum}' AFLAP_tmp/04/Call/${h}_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt)
			Cov=$(awk '{print $2}' AFLAP_tmp/04/Count/${h}_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt | sort -n | uniq -c | tail -n+2 | sort -nrk1,1 | head -n 1 | awk '{print $2}')
			if [[ $Cov == 1 ]]
				then
				Not1=0
				Is1=0
				Not1=$(awk '{print $2}' AFLAP_tmp/04/Count/${h}_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt | sort -n | uniq -c | awk 'NR > 5 {sum += $1}END{print sum}')
				Is1=$(awk '{print $2}' AFLAP_tmp/04/Count/${h}_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt | sort -n | uniq -c | awk 'NR == 2 {sum += $1}END{print sum}')
				if (( $Not1 > $Is1 ))
					then
					Cov=$(awk '{print $2}' AFLAP_tmp/04/Count/${h}_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt | sort -n | uniq -c | tail -n+6 | sort -nrk1,1 | head -n 1 | awk '{print $2}')
					else
					Cov=$(awk '{print $2}' AFLAP_tmp/04/Count/${h}_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt | sort -n | uniq -c | tail -n+2 | sort -nrk1,1 | head -n 1 | awk '{print $2}')
				fi
			fi
			echo -e "${h}\t${Mcou}\t${Cov}"
			else
			echo "Counts could not be found. Please rerun AFLAP"
			fi
		done | awk -v OFS='\t' '{if ($2 == 0) $3 = 0; print $0}' > AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerCount.txt
		Rscript $DIR/KmerCovXMarkerCount.R AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerCount.txt AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_KmerCovXMarkerCount.png
		LowCov=$(awk -v CC=$CovCut '$3 < CC' AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerCount.txt | wc -l)
		if (( $LowCov >= 1))
		then
		awk '$3 < 2 {print $1" appears to be low coverage, will be excluded"}' AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerCount.txt
#Loop through the ones which are high coverage and link files into a new directory.
#Build filtered genotype table
#Create header.
		fi
	cd  AFLAP_tmp/05/FilteredCall
	for v in `awk -v CC=$CovCut '$3 >= CC {print $1}' ../../../AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerCount.txt` 
	do
		ln -s ../../04/Call/${v}_${g}_m${mer}_L${Lo}_U${Up}_${P0}.txt .
	done
	cd ../
	awk -v OFS='\t' '{print $2, $1}' ../04/${g}_m${mer}_L${Lo}_U${Up}_$P0.Genotypes.MarkerID.tsv | paste - FilteredCall/* | awk -v SDU=$SegDistU -v SDL=$SegDistL '{for (i=3; i<=NF;i++) j+=$i; if(j/(NF-2) >= SDL && j/(NF-2) <= SDU) print $0 ; j=0 }' - > ${g}_m${mer}_L${Lo}_U${Up}_$P0.Filtered.Genotypes.MarkerID.tsv
	if (( $CovCut == 1 )); then 
	echo -e "AFLAP ran in low coverage mode. It is possible that two peaks will be shown in AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerSeg.png\nIf that is the case please rerun AFLAP.sh providing -d and -D for lower and upper limits for marker filtering."
	fi
	ls FilteredCall/* | sed "s/_${g}.*//" | sed 's/.*\///' > ${g}_m${mer}_L${Lo}_U${Up}_$P0.ProgHeader.txt
	cat <(cat ${g}_m${mer}_L${Lo}_U${Up}_$P0.ProgHeader.txt | tr '\n' '\t' | sed 's/\t$//' | awk -v OFS='\t' '{print "MarkerID", "MakerSeq", $0}') ${g}_m${mer}_L${Lo}_U${Up}_$P0.Filtered.Genotypes.MarkerID.tsv > ../../AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_$P0.GT.tsv
#Removing means isolates can be excluded by editing the Pedigree file
	rm FilteredCall/*
	cd ../../
	else
		echo "Draft genotype table for $g not found. Please rerun the pipeline"
		exit 1
	fi
done

exit
