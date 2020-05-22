#!/bin/bash -l

#Usage Statement

while getopts ':hM:1:2:W:p:' option; do
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

echo -e "\nBeginning AFLAP script 3/5"

#Check for infor from previous scripts.
if [[ -e AFLAP_tmp/01/LA.txt ]]
then
echo -e "\nLinkage analysis to be run using markers derived from:"
cat AFLAP_tmp/01/LA.txt
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

for g in 

#Could run below in a sub script to run both in parallel.
#./GetMarkers.sh -M $Mark1 -1 $Pare1 -2 $Pare2 -p $Out
#Would work well for providing the second with option W

#1. Identify unique mers to P1.
jellyfish query -s $Mark1 $Pare2 | awk '$2 == 0 {print ">"++i,"\n"$1}' > ./$t/${Out}Mark_P2zero.fa
#First instance of STDERR will overwrite file if present.
M1cou=$(grep -c '^>' ./$t/${Out}Mark_P2zero.fa)
echo "$M1cou mers unique to $Pare1 when compared to $Pare2" >&2

#2. Assemble unique
ABYSS -k25 -c0 -e0 ./$t/${Out}Mark_P2zero.fa -o ./$t/${Out}Mark_ABYSS.fa
MLcou=$(cat ./$t/${Out}Mark_ABYSS.fa | paste - - | awk '$2 < 61' | wc -l)
MLlen=$(cat ./$t/${Out}Mark_ABYSS.fa | paste - - | awk '$2 < 61{sum += $2}END{print sum}')
MMcou=$(cat ./$t/${Out}Mark_ABYSS.fa | paste - - | awk '$2 == 61' | wc -l)
MMlen=$(cat ./$t/${Out}Mark_ABYSS.fa | paste - - | awk '$2 == 61{sum += $2}END{print sum}')
MHcou=$(cat ./$t/${Out}Mark_ABYSS.fa | paste - - | awk '$2 > 61' | wc -l)
MHlen=$(cat ./$t/${Out}Mark_ABYSS.fa | paste - - | awk '$2 > 61{sum += $2}END{print sum}')
echo -e "$MLcou fragments smaller than 61bp, totalling $MLlen bp\n$MMcou fragments equal to 61bp, totalling $MMlen bp\n$MHcou fragments larger than 61bp, totalling $MHlen bp\n" >&2

#3. Filter Markers against hashes.
cat ./$t/${Out}Mark_ABYSS.fa | paste - - | awk '$2 == 61 {print $1"\n"substr($4,10,31)}' > ./$t/${Out}_61bp.fa
cat ./$t/${Out}Mark_ABYSS.fa | paste - - | awk '$2 > 61 {print $1"\n"substr($4,10,31)}' > ./$t/${Out}_62bp.fa
#Check P1 hash to ensure that it is between the boundaries.
#Boundaries can be parsed from header.
L=$(echo $Mark1 | sed 's/.*L//' | sed 's/_.*//')
U=$(echo $Mark1 | sed 's/.*U//' | sed 's/\..*//')
echo "Lower bound of $L obtained from file name" >&2
echo "Upper bound of $U obtained from file name" >&2
jellyfish query -s ./$t/${Out}_61bp.fa $Pare1 | awk -v Low="$L" -v Up="$U" '$2 >= Low && $2 <= Up {print ">"++i"\n"$1}' > ./$t/${Out}_61bp_het.fa
jellyfish query -s ./$t/${Out}_61bp_het.fa $Pare2 | awk '$2 == 0 {print ">"++i"\n"$1}' > ./$t/${Out}_61bp_MARKER.fa
jellyfish query -s ./$t/${Out}_62bp.fa $Pare1 | awk -v Low="$L" -v Up="$U" '$2 >= Low && $2 <= Up {print ">"++i"\n"$1}' > ./$t/${Out}_62bp_het.fa
jellyfish query -s ./$t/${Out}_62bp_het.fa $Pare2 | awk '$2 == 0 {print ">"++i"\n"$1}' > ./$t/${Out}_62bp_MARKER.fa

#4 Obtain ABYSS headers (Not currently essential)
#Obtain RevCom of ABYSS
cat ./$t/${Out}_61bp.fa | paste - - | cut -f 2 | tr 'ATCG' 'TAGC' | rev | paste <(cat ./$t/${Out}_61bp.fa | paste - - | cut -f 1) - | tr '\t' '\n' | cat - ./$t/${Out}_61bp.fa > ./$t/${Out}_61bp_RC.fa
cat ./$t/${Out}_61bp_MARKER.fa | paste - - | sort -k2,2 | join -1 2 -2 2 -o "1.1,0" - <(cat ./$t/${Out}_61bp_RC.fa | paste - - | sort -k2,2) | tr ' ' '\n' > ./${Out}_61bp_MARKER.fa
cat ./$t/${Out}_62bp.fa | paste - - | cut -f 2 | tr 'ATCG' 'TAGC' | rev | paste <(cat ./$t/${Out}_62bp.fa | paste - - | cut -f 1) - | tr '\t' '\n' | cat - ./$t/${Out}_62bp.fa > ./$t/${Out}_62bp_RC.fa
cat ./$t/${Out}_62bp_MARKER.fa | paste - - | sort -k2,2 | join -1 2 -2 2 -o "1.1,0" - <(cat ./$t/${Out}_62bp_RC.fa | paste - - | sort -k2,2) | tr ' ' '\n' > ./${Out}_62bp_MARKER.fa
M61=$(grep -c '^>' ./${Out}_61bp_MARKER.fa)
M62=$(grep -c '^>' ./${Out}_62bp_MARKER.fa)
echo -e "$M61 unique single copy markers equal to 61bp obtained from $Pare1 for downstream analysis" >&2
echo -e "$M62 unique single copy markers larger than 61 bp obtained from $Pare1 for downstream analysis" >&2

exit
