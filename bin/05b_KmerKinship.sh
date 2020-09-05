#!/bin/bash -l
#Obtain Kinship estimates. Further filter manually
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
#1 Check for Genotype table
#Again need a low vs high coverage flag
        P0=$(cat AFLAP_tmp/03/${g}_CrossedTo.txt | tr '\n' '_' | sed 's/_$//')
        Lo=$(awk -v var="$g" '$1 == var {print $2}' AFLAP_tmp/02/Boundaries.txt)
        Up=$(awk -v var="$g" '$1 == var {print $3}' AFLAP_tmp/02/Boundaries.txt)
	Rscript $DIR/KmerKinship.R AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_$P0.Filtered.Genotypes.MarkerID.tsv AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_$P0.ProgHeader.txt AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_$P0.KmerKinship.png AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_$P0.KmerKinship.txt
done
