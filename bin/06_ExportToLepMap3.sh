#!/bin/bash -l
#################################################
#       A shell script to export the genotype table to LepMap3
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

mkdir -p AFLAP_tmp/06
DIR=$(dirname $0)

#Inspect pedigeee file, decide if F1 or F2.
#Loop through per parent builing genotype file, LepMap3 format.
###########################################################
#	Notes on LepMap3 format
#	The format is a 10-field genotype, such that each allelic combination is specified by a space seperated number.
#	Decimals can be provided for likelihoods, here we are just working with absolutes.
###########################################################
#Provide additional header. 

#Loop through per parent.
for g in `cat AFLAP_tmp/01/LA.txt`
do
#Obtain coverage to get intermediate file labels
	if [[ -e AFLAP_tmp/03/${g}_CrossedTo.txt ]]
		then
		P0=$(cat AFLAP_tmp/03/${g}_CrossedTo.txt | tr '\n' '_' | sed 's/_$//')
		Lo=$(awk -v var="$g" '$1 == var {print $2}' AFLAP_tmp/02/Boundaries.txt)
		Up=$(awk -v var="$g" '$1 == var {print $3}' AFLAP_tmp/02/Boundaries.txt)
		else
		echo -e "Intermediate file is missing, please rerun AFLAP (s0603P1CT)"
		exit 1
	fi
	if [[ -e AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_${P0}.Filtered.Genotypes.MarkerID.tsv && -e AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_${P0}.ProgHeader.txt ]]
		then
#If F1 or F2
#Already validated that the file only contains one type of cross, so ok to look for first instance of parent and obtain cross information
		F=$(awk -v var=$g '$4 == var || $5 == var {print $2}' $Ped | head -n 1)
		echo -e "$Ped defines the cross as F$F"
		Z=$(awk -v var=$g '$4 == var || $5 == var {if ($4 == var) print 4; else if ($5 == var) print 5}' $Ped | head -n 1)
#Check to make sure P1 isn't in both columns already performed?
		if (( $Z == 3 ))
			then
			echo -e "Uh oh, $g doesn't seem to be correctly specified in $Ped. AFLAP should not have advanced this far. Please validate $Ped and try rerunning AFLAP (S0600PPPf)"
			exit 1
			else
			echo -e "$g detected in column $Z of $Ped, this will be reflected in the output."
		fi
#If parent in column one or column 2
###########################################################
#	Following will loop through the header file generated in script 5 and output in a LepMap3 format
#	Will modify the genotype scores to be suitable for LepMap3
###########################################################

		if (( $Z == 4 && $F == 1 ))
			then
			cat \
			<(cat \
			<(echo -e "CHR CHR CHR CHR CHR CHR\nPOS POS POS POS POS POS\n${g}x${P0} $g 0 0 1 0\n${g}x${P0} $P0 0 0 2 0") \
			<(awk -v P1=$g -v P2=$P0 '{print P1"x"P2, $1, P1, P2, 0, 0}' AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_${P0}.ProgHeader.txt) \
			| awk -f $DIR/Transpose.awk -) \
			<(paste \
			<(awk -v OFS='\t' '{print $1, $2}'  AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_${P0}.Filtered.Genotypes.MarkerID.tsv) \
			<(cut -f 3- AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_${P0}.Filtered.Genotypes.MarkerID.tsv | awk -v OFS='\t' '{print 1,0,$0}' | sed 's/$/\t/' | sed 's/1\t/0 1 0 0 0 0 0 0 0 0|\t/g' | 
			sed 's/0\t/1 0 0 0 0 0 0 0 0 0|\t/g' | sed 's/|//g' | sed 's/\t$//')\
			)> AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}.ForLepMap3.tsv
		elif (( $Z == 5 && $F == 1 ))
			then
			cat \
			<(cat \
			<(echo -e "CHR CHR CHR CHR CHR CHR\nPOS POS POS POS POS POS\n${P0}x${g} $P0 0 0 1 0\n${P0}x${g} $g 0 0 2 0") \
			<(awk -v P1=$P0 -v P2=$g '{print P1"x"P2, $1, P1, P2, 0, 0}' AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_${P0}.ProgHeader.txt) \
			| awk -f $DIR/Transpose.awk -) \
			<(paste \
			<(awk -v OFS='\t' '{print $1, $2}'  AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_${P0}.Filtered.Genotypes.MarkerID.tsv) \
			<(cut -f 3- AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_${P0}.Filtered.Genotypes.MarkerID.tsv | awk -v OFS='\t' '{print 0,1,$0}' | sed 's/$/\t/' | sed 's/1\t/0 1 0 0 0 0 0 0 0 0|\t/g' | 
			sed 's/0\t/1 0 0 0 0 0 0 0 0 0|\t/g' | sed 's/|//g' | sed 's/\t$//')\
			)> AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}.ForLepMap3.tsv
		elif (( $Z == 4 && $F == 2 ))
			then
			cat \
			<(cat \
			<(echo -e "CHR CHR CHR CHR CHR CHR\nPOS POS POS POS POS POS\n${g}x${P0} $g 0 0 1 0\n${g}x${P0} $P0 0 0 2 0\n${g}x${P0} DUM1 $g $P0 1 0\n${g}x${P0} DUM2 $g $P0 2 0") \
			<(awk '{print P1"x"P2, $1, "DUM1", "DUM2", 0, 0}' AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_${P0}.ProgHeader.txt) \
			| awk -f $DIR/Transpose.awk -) \
			<(paste \
			<(awk -v OFS='\t' '{print $1, $2}'  AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_${P0}.Filtered.Genotypes.MarkerID.tsv) \
			<(cut -f 3- AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_${P0}.Filtered.Genotypes.MarkerID.tsv | awk -v OFS='\t' '{print 1,0,2,2,$0}' | sed 's/$/\t/' | sed 's/1\t/0 0 0 0 1 0 0 0 0 0|\t/g' | 
			sed 's/0\t/1 0 0 0 0 0 0 0 0 0|\t/g' | sed 's/2\t/0 1 0 0 0 0 0 0 0 0|\t/g' | sed 's/|//g' | sed 's/\t$//')\
			)> AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}.ForLepMap3.tsv
		elif (( $Z == 5 && $F == 2 ))
			then
			cat \
			<(cat \
			<(echo -e "CHR CHR CHR CHR CHR CHR\nPOS POS POS POS POS POS\n${P0}x${g} $P0 0 0 1 0\n${P0}x${g} $g 0 0 2 0\n${g}x${P0} DUM1 $P0 $g 1 0\n${g}x${P0} DUM2 $P0 $g 2 0") \
			<(awk '{print P1"x"P2, $1, "DUM1", "DUM2", 0, 0}' AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_${P0}.ProgHeader.txt) \
			| awk -f $DIR/Transpose.awk -) \
			<(paste \
			<(awk -v OFS='\t' '{print $1, $2}'  AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_${P0}.Filtered.Genotypes.MarkerID.tsv) \
			<(cut -f 3- AFLAP_tmp/05/${g}_m${mer}_L${Lo}_U${Up}_${P0}.Filtered.Genotypes.MarkerID.tsv | awk -v OFS='\t' '{print 0,1,2,2,$0}' | sed 's/$/\t/' | sed 's/1\t/0 0 0 0 1 0 0 0 0 0|\t/g' | 
			sed 's/0\t/1 0 0 0 0 0 0 0 0 0|\t/g' | sed 's/2\t/0 1 0 0 0 0 0 0 0 0|\t/g' | sed 's/|//g' | sed 's/\t$//')\
			)> AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}.ForLepMap3.tsv
		else
			echo -e "Oh no! $Ped doesn't seem to be formatted correctly. AFLAP should not have advanced this far. Please validate $Ped and try rerunning AFLAP (S0600Pf)"
		fi
		else
		echo -e "Intermediate file is missing, please rerun AFLAP (s0605FGTT)"
		exit 1
	fi
#Print location of LepMap3 compatible genotype table
PWD=$(pwd)
echo -e "LepMap3 Genotype table for ${g}:\n$PWD/AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}.ForLepMap3.tsv"
done

exit 

