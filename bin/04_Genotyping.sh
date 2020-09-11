#!/bin/bash -l
#################################################
#	A shell script to call genotypes of progeny using makrers derived from a parent and progeny JELLYFISH hashes
#	For optimal calls, a k-mer should be observed twice.
#################################################
Name=$(basename $0)
usage="${Name}; [-h] [-P] [-m] -- A script to genotype progeny (AFLAP 4/6).
Options
	-h show this help message
	-P Pedigree file, required. See AFLAP README for more information.
	-m K-mer size. Optional. Default [31]
Temporary files will be output to AFLAP_tmp/04."
#Option block
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

#0. Dependency Check
#a JELLYFISH
JF=$(jellyfish --version)
if [[ $JF =~ ^jellyfish ]]; then echo "$JF detected"; else echo "jellyfish not detected, please modify your PATH" ; exit 1 ; fi

mkdir -p AFLAP_tmp/04/Count
mkdir -p AFLAP_tmp/04/Call

#Check For markers
for g in `cat AFLAP_tmp/01/LA.txt`
do
	echo -e "Calling GT for $g derived markers"
	#Identify correct file based on intermediate results.
	P0=$(cat AFLAP_tmp/03/${g}_CrossedTo.txt | tr '\n' '_' | sed 's/_$//')
	Lo=$(awk -v var="$g" '$1 == var {print $2}' AFLAP_tmp/02/Boundaries.txt)
	Up=$(awk -v var="$g" '$1 == var {print $3}' AFLAP_tmp/02/Boundaries.txt)
	if [[ -e AFLAP_tmp/03/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa ]]
	then 
		Mcou=$(grep -c '^>' AFLAP_tmp/03/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa)
		echo -e "$Mcou markers identified in AFLAP_tmp/03/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa. These will be surveyed against progeny"
	else
		echo -e "Can't automatically locate marker file. Please try rerunning AFLAP.sh."
	fi
	for h in `awk -v var="$g" '$4 == var || $5 == var {print $1}' $Ped | sort -u`
	do
		if [[ -e AFLAP_tmp/01/ProgCounts/$h.jf${mer} ]]
		then
			if [[ -e AFLAP_tmp/04/Count/${h}_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt && -e AFLAP_tmp/04/Call/${h}_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt ]]; then
			echo -e "Genotype of $g markers for $h previously calulated, will use these"
			else
			jellyfish query -s AFLAP_tmp/03/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa AFLAP_tmp/01/ProgCounts/$h.jf${mer} > AFLAP_tmp/04/Count/${h}_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt
#Need a low cov and high cov option
			awk '{if($2 > 1) print 1; else print 0}' AFLAP_tmp/04/Count/${h}_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt > AFLAP_tmp/04/Call/${h}_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt
			fi
		else
			echo -e "No hash for $h detected. Please rerun 01_JELLYFISH.sh"
			exit 1
		fi
	done
	echo -e "GT calling for $g derived markers complete"
	awk '{print $1}' AFLAP_tmp/04/Count/${h}_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt | paste - AFLAP_tmp/04/Call/*_${g}_m${mer}_L${Lo}_U${Up}_$P0.txt > AFLAP_tmp/04/${g}_m${mer}_L${Lo}_U${Up}_$P0.Genotypes.tsv
	cat AFLAP_tmp/03/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa | paste - - | sed 's/>//' | sort -k2,2 | join -1 2 - <(sort -k1,1 AFLAP_tmp/04/${g}_m${mer}_L${Lo}_U${Up}_$P0.Genotypes.tsv) | tr ' ' '\t' > AFLAP_tmp/04/${g}_m${mer}_L${Lo}_U${Up}_$P0.Genotypes.MarkerID.tsv
	rm AFLAP_tmp/04/${g}_m${mer}_L${Lo}_U${Up}_$P0.Genotypes.tsv
done
exit
