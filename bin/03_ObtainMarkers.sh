#!/bin/bash -l
#################################################
#	A shell script to derive single copy k-mers that are unique to a parent. These are then used a markers.
#	To reduce redundancy k-mers are assembled using ABySS.
#	To enable the use of a consistent has size, the markers are reduced to a sequnce length equal to option m.
#	This consistent marker length then only need to be surveyed against only one progeny hash.
#################################################
Name=$(basename $0)
usage="${Name}; [-h] [-P] [-m] -- A script to obtain single copy k-mers from parental JELLYFISH hashes (AFLAP 3/6).
Options
	-h show this help message
	-P Pedigree file, required. See AFLAP README for more information.
	-m K-mer size. Optional. Default [31]
Temporary files will be output to AFLAP_tmp/03."
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
mkdir -p AFLAP_tmp/03/ParentalMarkers

#Strip '#' from pedigree file
awk '$0 !~ /#/' $Ped > AFLAP_tmp/Pedigree.txt
Ped=AFLAP_tmp/Pedigree.txt

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
		if [[ -e AFLAP_tmp/02/ParentalHisto/${g}_m${mer}_L${Lo}_U${Up}.fa ]]
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
	cp AFLAP_tmp/02/ParentalHisto/${g}_m${mer}_L${Lo}_U${Up}.fa AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}.fa
	#Copied into 03 so we can overwrite in the next stage.
#3. Filter against opposing parents (another loop?)
	P0=$(cat AFLAP_tmp/03/${g}_CrossedTo.txt | tr '\n' '_' | sed 's/_$//')
	if [[ -e AFLAP_tmp/03/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa ]]
		then
		echo -e "Previously calculated markers detected. Want to calculate new markers, please deleted:\n./AFLAP_tmp/03/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa\nSummary of markers available:\n./$g.MarkerReport.txt"
		else
	for f in `cat AFLAP_tmp/03/${g}_CrossedTo.txt`
	do
	# Identify opposing parental hashes.
		if  [[ -e AFLAP_tmp/01/ParentalCounts/$f.jf${mer} ]]
		then
		echo "Intersecting $g with $f"
	# Filter and overwrite input.
		jellyfish query -s AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}.fa AFLAP_tmp/01/ParentalCounts/$f.jf${mer} | awk '$2 == 0 {print ">"++i,"\n"$1}' > AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}.txt
		mv AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}.txt AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}.fa
		Mlcou=$(grep -c '^>' AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}.fa)
		echo "${Mlcou} $g ${mer}-mers remain after filtering against $f"
		else
		echo -e "JELLYFISH results of $f can't be located. Please rerun 01_JELLYFISH.sh."
		exit 1
		fi
	done
#4. Assemble
#Optimal k for 31-mer=25. For 21-mer=19
	if (( $mer == 31)); then k=25 ; elif (($mer == 21)); then k=19 ; else k=$(echo $mer-2 | bc) ; fi 
	echo -e "Running ABySS with -k set to $k"
	ABYSS -k$k -c0 -e0 AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}.fa -o AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark1.fa
#5. Extract
	ak=$(echo ${mer}+${mer}-1 | bc)
	cat AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark1.fa | paste - - | awk -v mer=$mer -v ak=$ak '$2 >= ak {print $1"_"$2"\n"substr($4,10,mer)}' > AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark2.fa
#6. Refilter against self for mers between boundaries.
	jellyfish query -s AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark2.fa AFLAP_tmp/01/ParentalCounts/$g.jf${mer} | awk -v Low="$Lo" -v Up="$Up" '$2 >= Low && $2 <= Up {print ">"++i"\n"$1}' > AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark3.fa
#7. Refilter against other parents (another loop).
	for f in `cat AFLAP_tmp/03/${g}_CrossedTo.txt`
	do
		jellyfish query -s AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark3.fa AFLAP_tmp/01/ParentalCounts/$f.jf${mer} | awk '$2 == 0 {print ">"++i"\n"$1}' > AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark4.fa
		mv AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark4.fa AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark3.fa
	done
#8. Export final marker set with ABySS conserved headers.
	cat AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark2.fa | paste - - | cut -f 2 | tr 'ATCG' 'TAGC' | rev | paste <(cat AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark2.fa | paste - - | cut -f 1) - | tr '\t' '\n' | cat - AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark2.fa > AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark2RC.fa
	cat AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark3.fa | paste - - | sort -k2,2 | join -1 2 -2 2 -o "2.1,0" - <(cat AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark2RC.fa | paste - - | sort -k2,2) | tr ' ' '\n' > AFLAP_tmp/03/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa
#9. Export stats.
	FragCou=$(grep -c '^>' AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark1.fa)
	Cou61=$(cat AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark1.fa | paste - - | awk -v ak=$ak '$2 == ak' | wc -l)
	Cou62=$(cat AFLAP_tmp/03/${g}_m${mer}_L${Lo}_U${Up}_Mark1.fa | paste - - | awk -v ak=$ak '$2 > ak' | wc -l)
	MarCou=$(grep -c '^>' AFLAP_tmp/03/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa)
	Mar61=$(cat AFLAP_tmp/03/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa | paste - - | sed 's/_/ /' | awk -v ak=$ak '$2 == ak' | wc -l)
	Mar62=$(cat AFLAP_tmp/03/ParentalMarkers/${g}_m${mer}_MARKERS_L${Lo}_U${Up}_$P0.fa | paste - - | sed 's/_/ /' | awk -v ak=$ak '$2 > ak' | wc -l)
	echo -e "Report for $g
Number of $mer-mers input into assembly: $Mlcou
Number of fragments assembled: $FragCou
Number of fragments == $ak bp: $Cou61
Number of fragments > $ak bp: $Cou62
Number of markers after refiltering: $MarCou
Number of markers == $ak bp: $Mar61
Number of markers > $ak bp: $Mar62
" | tee $g.MarkerReport.txt
	fi
done

exit
