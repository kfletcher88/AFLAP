#!/bin/bash
Name=$(basename $0)
usage="${Name}; [-h] [-P] [-t] [-m] -- A script to run all stages of AFLAP.

Options
        -h show this help message
        -P Pedigree file, required. See AFLAP README for more information.
        -m K-mer size. Optional. Default [31]
        -t Threads for JELLYFISH counting. Optional. Default [4]
	-r Individual to remove. All other options will be ignored."

while getopts ':khP:t:m:r:L:' option; do
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
		k)  kin=1
			;;
		L)  LOD=$OPTARG
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
if (( $kin == 1 )); then echo -e "\n\e[31mRunning Kmer kinship\e[0m" ; $DIR/bin/05b_KmerKinship.sh -P $Ped -m $mer ; fi &&
echo -e "\n\e[31mBeginning Step 6/6\e[0m" &&
$DIR/bin/06_ExportToLepMap3.sh -P $Ped -m $mer &&
if [[ -z $LOD ]]; then echo -e "No LOD cutoffs provided, so AFLAP will not run LepMap3. LOD cutoffs can be provided with the -L flag" ;
else
for g in `echo $LOD | tr ',' '\n'` ; do $DIR/bin/07_LepMap3.sh -P $Ped -m $mer -T $thread -L $g ; done 
fi
exit
else
echo -e "Removing intermediate progeny files generated for $rem"
find AFLAP_Intermediate/ -name "${rem}*" -exec rm {} +
find AFLAP_tmp/ -name "${rem}*" -exec rm {} +
exit
fi
#Add additional scripts for running LepMap3
