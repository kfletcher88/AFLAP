#!/bin/bash -l
#Obtain Segregation Histograms

#!/bin/bash -l

#Usage Statement

while getopts ':hP:m:' option; do
        case "$option" in
                h)  echo "$usage"
                         exit
                         ;;
                P)  Ped=$OPTARG
                         ;;
                m)  mer=$OPTARG
                         ;;
                \?) printf "illegal option: -%s\n\n" "$OPTARG" >&2
                    echo "$usage"
                    exit 1
                         ;;
        esac
done

if [[ -e AFLAP_tmp/01/LA.txt && -e AFLAP_tmp/02/Boundaries.txt ]]
then
echo -e "\nIntermediate files detected"
else
echo -e "\n Could not find output of previous scripts. Please rerun full pipeline."
exit 1
fi

mkdir -p AFLAP_Intermediate/SegregationInformation
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
		ProC=$(awk 'NR == 1 {print NF-2}' AFLAP_tmp/04/${g}_m${mer}_L${Lo}_U${Up}_$P0.Genotypes.MarkerID.tsv)
		echo "Genotype table for $g detected. Summarizing"
		sed 's/_/ /' AFLAP_tmp/04/${g}_m${mer}_L${Lo}_U${Up}_$P0.Genotypes.MarkerID.tsv | awk -v var="$kk" '$3 == var {for (i=4; i<=NF;i++) j+=$i; print j; j=0 }' | sort -n | uniq -c | awk -v var=$ProC '{print $2/var, $1}' > AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerEqual${kk}.hist
		sed 's/_/ /' AFLAP_tmp/04/${g}_m${mer}_L${Lo}_U${Up}_$P0.Genotypes.MarkerID.tsv | awk -v var="$kk" '$3 > var {for (i=4; i<=NF;i++) j+=$i; print j; j=0 }' | sort -n | uniq -c | awk -v var=$ProC '{print $2/var, $1}' > AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerOver${kk}.hist
		sed 's/_/ /' AFLAP_tmp/04/${g}_m${mer}_L${Lo}_U${Up}_$P0.Genotypes.MarkerID.tsv | awk '{for (i=4; i<=NF;i++) j+=$i; print j; j=0 }' | sort -n | uniq -c | awk -v var=$ProC '{print $2/var, $1}' > AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_AllMarkers.hist
		Rscript $DIR/bin/SegStats.R AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerEqual${kk}.hist AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerOver${kk}.hist  AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_AllMarkers.hist AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerSeg.png
		for h in `awk -v var="$g" '$4 == var || $5 == var {print $1}' $Ped | sort -u`
		do
			if [[ -e AFLAP_tmp/04/Call/${h}_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt ]]
			then
			Mcou=$(awk '{sum += $1}END{print sum}' AFLAP_tmp/04/Call/${h}_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt)
			Cov=$(awk '{print $2}' AFLAP_tmp/04/Count/${h}_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt | sort -n | uniq -c | tail -n+2 | sort -nrk1,1 | head -n 1 | awk '{print $2}')
			if [[ $Cov == 1 ]]
				then
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
		Rscript $DIR/bin/KmerCovXMarkerCount.R AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerCount.txt AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_KmerCovXMarkerCount.png
#Something like this should be placed before the genotype table, so low coverage isolates are not incorporated.
		LowCov=$(awk '$3 < 2' AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerCount.txt | wc -l)
		if (( $LowCov >= 1))
		then
		awk '$3 < 3 {print $1" appears to be low coverage, will be excluded"}' AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}_MarkerCount.txt
#Loop through the ones which are high coverage and link files into a new directory.
#Build filtered genotype table
#Create header.
		fi
		
	else
		echo "Genotype table for $g not found. Please rerun the pipeline"
		exit 1
	fi
done

exit
