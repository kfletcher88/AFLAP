#!/bin/bash
Name=$(basename $0)
usage="${Name}; [-h] [-P] [-t] [-m] -- A script to run all stages of AFLAP.

Options
        -h show this help message
        -P Pedigree file, required. See AFLAP README for more information.
        -m K-mer size. Optional. Default [31]
        -t Threads for JELLYFISH counting. Optional. Default [4]
	-r Individual to remove. All other options will be ignored."

while getopts ':hP:t:m:r:' option; do
        case "$option" in
                h)  echo "$usage"
                         exit
                         ;;
                P)  Ped=$OPTARG
                         ;;
                t)  thread=$OPTARG
                        ;;
                m)  mer=$OPTARG
                        ;;
		r)  rem=$OPTARG
			;;
                \?) printf "illegal option: -%s\n\n" "$OPTARG" >&2
                    echo "$usage"
                    exit 1
                         ;;
        esac
done

if [[ -z $rem ]]; then
DIR=$(dirname $0)
#00 Dependency check
echo -e "\n\e[31mBeginning Step 1/6\e[0m" &&
$DIR/bin/01_JELLYFISH.sh -P $Ped -t $thread -m $mer &&
echo -e "\n\e[31mBeginning Step 2/6\e[0m" &&
$DIR/bin/02_ExtractSingleCopyMers.sh -P $Ped -m $mer &&
echo -e "\n\e[31mBeginning Step 3/6\e[0m" &&
$DIR/bin/03_ObtainMarkers.sh -P $Ped -m $mer &&
echo -e "\n\e[31mBeginning Step 4/6\e[0m" &&
$DIR/bin/04_Genotyping.sh -P $Ped -m $mer &&
echo -e "\n\e[31mBeginning Step 5/6\e[0m" &&
$DIR/bin/05_ObtainSegStats.sh -P $Ped -m $mer &&
echo -e "\n\e[31mBeginning Step 6/6\e[0m" &&
$DIR/bin/06_ExportToLepMap3.sh -P $Ped -m $mer &&
exit
else
echo -e "Removing intermediate progeny files generated for $rem"
find AFLAP_Intermediate/ -name "${rem}*" -exec rm {} +
find AFLAP_tmp/ -name "${rem}*" -exec rm {} +
exit
fi
#Add additional scripts for running LepMap3
