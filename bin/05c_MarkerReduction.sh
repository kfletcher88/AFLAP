#!/bin/bash -l
#################################################
#       Optional script to reduce the number of markers output in the final genotype table, reducing runtime.
#################################################
Name=$(basename $0)
usage="${Name}; [-h] [-P] [-m] [-U]-- A script to genotype progeny (AFLAP 5b/6).
Options
        -h show this help message
        -P Pedigree file, required. See AFLAP README for more information.
        -m K-mer size. Optional. Default [31]
	-U Maximum number of markers to output."
#Option block
while getopts 'P:m:U:' option; do
        case "$option" in
                P)  Ped=$OPTARG
                         ;;
                m)  mer=$OPTARG
                         ;;
		U)  Max=$OPTARG
			 ;;
                \?) exit 1
                         ;;
        esac
done
if [[ -z $Ped || -z $mer || -z $Max ]]; then echo "Weird, are you trying to run this apart from AFLAP.sh. That is not recommended. Please use AFLAP.sh -U to launch correctly"; exit 1 ; fi
if [[ -e AFLAP_tmp/01/LA.txt && -e AFLAP_tmp/02/Boundaries.txt ]]
then
echo -e "\nIntermediate files detected"
else
echo -e "\n Could not find output of previous scripts. Please rerun full pipeline."
exit 1
fi
#Find bin directory for Rscripts
DIR=$(dirname $0)

for g in `cat AFLAP_tmp/01/LA.txt`
do
#Check for Genotype table
        P0=$(cat AFLAP_tmp/03/${g}_CrossedTo.txt | tr '\n' '_' | sed 's/_$//')
        Lo=$(awk -v var="$g" '$1 == var {print $2}' AFLAP_tmp/02/Boundaries.txt)
        Up=$(awk -v var="$g" '$1 == var {print $3}' AFLAP_tmp/02/Boundaries.txt)
	MarkerCount=$(cat AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_$P0.Filtered.Genotypes.MarkerID.tsv | wc -l)
	if (( $MarkerCount < $Max )); then echo -e "\n$MarkerCount markers detetected in genotype file. This is less than $Max so no action taken"
	else mv AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_$P0.Filtered.Genotypes.MarkerID.tsv AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_$P0.Filtered.Genotypes.${MarkerCount}_MarkerID.tsv
	shuf -n $Max AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_$P0.Filtered.Genotypes.${MarkerCount}_MarkerID.tsv > AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_$P0.Filtered.Genotypes.MarkerID.tsv
	fi
done
