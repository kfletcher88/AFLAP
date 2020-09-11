#!/bin/bash -l
#################################################
#	Optional script to plot a heatmap, approximating kinship
#	The script will not exclude isolates from downstram analysis, but may guide the user in which progeny, if any to remove
#################################################
Name=$(basename $0)
usage="${Name}; [-h] [-P] [-m] -- A script to genotype progeny (AFLAP 5b/6).
Options
        -h show this help message
        -P Pedigree file, required. See AFLAP README for more information.
        -m K-mer size. Optional. Default [31]
Temporary files will be output to AFLAP_tmp/KmerKin."
#Option block
while getopts 'P:m:' option; do
        case "$option" in
                P)  Ped=$OPTARG
                         ;;
                m)  mer=$OPTARG
                         ;;
                \?) exit 1
                         ;;
        esac
done
if [[ -z $Ped || -z $mer ]]; then echo "Weird, are you trying to run this apart from AFLAP.sh. That is not recommended. Please use AFLAP.sh -k to launch correctly"; exit 1 ; fi
if [[ -e AFLAP_tmp/01/LA.txt && -e AFLAP_tmp/02/Boundaries.txt ]]
then
echo -e "\nIntermediate files detected"
else
echo -e "\n Could not find output of previous scripts. Please rerun full pipeline."
exit 1
fi

mkdir -p AFLAP_tmp/KmerKin
#Find bin directory for Rscripts
DIR=$(dirname $0)

for g in `cat AFLAP_tmp/01/LA.txt`
do
#Check for Genotype table
        P0=$(cat AFLAP_tmp/03/${g}_CrossedTo.txt | tr '\n' '_' | sed 's/_$//')
        Lo=$(awk -v var="$g" '$1 == var {print $2}' AFLAP_tmp/02/Boundaries.txt)
        Up=$(awk -v var="$g" '$1 == var {print $3}' AFLAP_tmp/02/Boundaries.txt)
	Rscript $DIR/KmerKinship.R AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_$P0.Filtered.Genotypes.MarkerID.tsv AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_$P0.ProgHeader.txt AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_$P0.KmerKinship.png AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_$P0.KmerKinship.txt
done
