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

t=AFLAP_tmp
mkdir -p $t

#Check arguments.
if [[ -z $Ped ]]; then
        echo -e "ERROR: Pedigree file not provided\n"
        echo "$usage"
        exit 1
fi
if [[ -z $mer ]]; then
echo "mer size not specified, will proceed with default [31]"
mer=31
fi

#Check for required tmp files from previous scripts.
if [[ -e AFLAP_tmp/01/LA.txt ]]
then
echo -e "\nMarkers to be derived from:"
cat AFLAP_tmp/01/LA.txt
else
echo -e "\n Could not find output of previous scripts. Please rerun 01_JELLYFISH.sh."
exit 1
fi

if [[ -e AFLAP_tmp/01/Crosses.txt ]]
then
echo -e "Crosses to be analyzed include:"
awk '{print $3"x"$4}' AFLAP_tmp/01/Crosses.txt | head
else
echo -e "\n Could not find output of previous scripts. Please rerun 01_JELLYFISH.sh."
exit 1
fi

if [[ -e AFLAP_tmp/02/Boundaries.txt ]]
then
echo -e "${mer}-mer boundaries identified from previous run:"
head AFLAP_tmp/02/Boundaries.txt
else
echo -e "\n Could not find output of previous scripts. Please rerun 02_ExtractSingleCopyMers.sh."
exit 1
fi

#Make new tmp and Intermdiate directories:
mkdir -p AFLAP_tmp/03
mkdir -p AFLAP_Intermediate/ParentalMarkers

#0. Dependency Check

#a. ABYSS
AB=$(ABYSS --version | head -n 1)
if [[ $AB =~ ^ABYSS ]]; then echo "$AB detected"; else echo "ABySS not detected, please modify your PATH" ; exit 1 ; fi

#b JELLYFISH
JF=$(jellyfish --version)
if [[ $JF =~ ^jellyfish ]]; then echo "$JF detected"; else echo "jellyfish not detected, please modify your PATH" ; exit 1 ; fi

#Only assemble markers for those which boundaries were identified.
for g in `cat AFLAP_tmp/01/LA.txt`
do
	echo -e "Begining analysis for $g"
#1. Identify upper and lower boundaries for file ID and filtering.
	echo -e "Identifying ${mer}-mer boundaries used previously."
	Lo=$(awk -v var="$g" '$1 == var {print $2}' AFLAP_tmp/02/Boundaries.txt)
	Up=$(awk -v var="$g" '$1 == var {print $3}' AFLAP_tmp/02/Boundaries.txt)
	echo -e "Lower boundary set to $Lo"
	echo -e "Upper boundary set to $Up"
		if [[ -e AFLAP_Intermediate/ParentalHisto/${g}_m${mer}_L${Lo}_U${Up}.fa ]]
		then
		echo -e "Redults with these boundaries detected. Good to proceed."
		else
		echo -e "Hmm can't find results with these boundaries. Please try to rerun the pipeline."
		exit 1
		fi
#2. Identify opposing parental hashes to filter against.
	echo -e "Identifyng parents crossed to $g"
	awk '{print $3, $4}' AFLAP_tmp/01/Crosses.txt | awk -v var="$g" '$0 ~ var {print $1"\n"$2}' | awk -v var="$g" '$1 != var' > AFLAP_tmp/03/${g}_CrossedTo.txt
	echo -e "$g identified as crossed to:"
	head AFLAP_tmp/03/${g}_CrossedTo.txt
	cp AFLAP_Intermediate/ParentalHisto/${g}_m${mer}_L${Lo}_U${Up}.fa AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}.fa
	#Copied into tmp so we can overwrite in the next stage.
#3. Filter against opposing parents (another loop?)
	P0=$(cat AFLAP_tmp/03/${g}_CrossedTo.txt | tr '\n' '_' | sed 's/_$//')
	if [[ -e AFLAP_Intermediate/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa ]]
		then
		echo -e "Previously calculated markers detected. Want to calculate new markers, please deleted:\n./AFLAP_Intermediate/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa\nSummary of markers available:\n./$g.MarkerReport.txt"
		else
	for f in `cat AFLAP_tmp/03/${g}_CrossedTo.txt`
	do
	# Identify opposing parental hashes.
		if  [[ -e AFLAP_Intermediate/ParentalCounts/$f.jf${mer} ]]
		then
		echo "Intersecting $g with $f"
	# Filter and overwrite input.
		jellyfish query -s AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}.fa AFLAP_Intermediate/ParentalCounts/$f.jf${mer} | awk '$2 == 0 {print ">"++i,"\n"$1}' > AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}.txt
		mv AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}.txt AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}.fa
		Mlcou=$(grep -c '^>' AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}.fa)
		echo "${Mlcou} $g ${mer}-mers remain after filtering against $f"
		else
		echo -e "JELLYFISH results of $f can't be located. Please rerun 01_JELLYFISH.sh."
		exit 1
		fi
	done
#4. Assemble
	ABYSS -k25 -c0 -e0 AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}.fa -o AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark1.fa
#5. Extract
	cat AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark1.fa | paste - - | awk '$2 >= 61 {print $1"_"$2"\n"substr($4,10,31)}' > AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark2.fa
#6. Refilter against self for mers between boundaries.
	jellyfish query -s AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark2.fa AFLAP_Intermediate/ParentalCounts/$g.jf${mer} | awk -v Low="$Lo" -v Up="$Up" '$2 >= Low && $2 <= Up {print ">"++i"\n"$1}' > AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark3.fa
#7. Refilter against other parents (another loop).
	for f in `cat AFLAP_tmp/03/${g}_CrossedTo.txt`
	do
		jellyfish query -s AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark3.fa AFLAP_Intermediate/ParentalCounts/$f.jf${mer} | awk '$2 == 0 {print ">"++i"\n"$1}' > AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark4.fa
		mv AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark4.fa AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark3.fa
	done
#8. Export final marker set with ABySS conserved headers.
	cat AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark2.fa | paste - - | cut -f 2 | tr 'ATCG' 'TAGC' | rev | paste <(cat AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark2.fa | paste - - | cut -f 1) - | tr '\t' '\n' | cat - AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark2.fa > AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark2RC.fa
	cat AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark3.fa | paste - - | sort -k2,2 | join -1 2 -2 2 -o "2.1,0" - <(cat AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark2RC.fa | paste - - | sort -k2,2) | tr ' ' '\n' > AFLAP_Intermediate/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa
#9. Export stats.
	FragCou=$(grep -c '^>' AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark1.fa)
	Cou61=$(cat AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark1.fa | paste - - | awk '$2 == 61' | wc -l)
	Cou62=$(cat AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark1.fa | paste - - | awk '$2 > 61' | wc -l)
	MarCou=$(grep -c '^>' AFLAP_Intermediate/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa)
	Mar61=$(cat AFLAP_Intermediate/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa | paste - - | sed 's/_/ /' | awk '$2 == 61' | wc -l)
	Mar62=$(cat AFLAP_Intermediate/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa | paste - - | sed 's/_/ /' | awk '$2 > 61' | wc -l)
	echo -e "Report for $g
Number of $mer-mers input into assembly: $Mlcou
Number of fragments assembled: $FragCou
Number of fragments == 61 bp: $Cou61
Number of fragments > 61 bp: $Cou62
Number of markers after refiltering: $MarCou
Number of markers == 61 bp: $Mar61
Number of markers > 61 bp: $Mar62
" | tee $g.MarkerReport.txt
	fi
done

exit
