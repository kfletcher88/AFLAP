#!/bin/bash
Name=$(basename $0)
usage="${Name}; [-h] [-P] [-t] [-m] -- A script to run all stages of AFLAP.

Options
        -h show this help message
        -P Pedigree file, required. See AFLAP README for more information.
        -m K-mer size. Optional. Default [31]
        -t Threads for JELLYFISH counting. Optional. Default [4]
	-r Individual to remove. All other options will be ignored.
	-L LOD score - Will run LepMap3 with minimum LOD.
	-d Lower boundary for marker cut off. Can be used to filter for segregation distortion [0.2].
	-D Upper boundary for marker cut off. Can be used to filter for segregation distortion [0.8].
	-k Run kinship estimation.
	-x Run with low coverage parameters.
	-U Maximum number of markers to output in the genotype tables output under ./AFLAP_Results/"

while getopts ':kxhP:t:m:r:L:d:D:U:' option; do
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
		x)  LowCov=1
			;;
		d)  Sdl=$OPTARG
			;;
		D)  sdu=$OPTARG
			;;
		U)  Max=$OPTARG
			;;
                \?) printf "illegal option: -%s\n\n" "$OPTARG" >&2
                    echo "$usage"
                    exit 1
                         ;;
        esac
done

if [[ -z $rem && -z $Ped ]]; then echo -e "Either pedigree file (-P), or progeny individual intermediate data to remove (-r) required" ; echo "$usage" ; exit 1 ; fi
if [[ -z $thread ]]; then
echo "Threads not specified, will proceed with default [4]"
thread=4
fi
if [[ -z $Sdl ]]; then Sdl=0.2 ; fi
if [[ -z $Sdu ]]; then Sdu=0.8 ; fi
if [[ -z $mer ]]; then
echo "mer size not specified, will proceed with default [31]"
mer=31
fi
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
if [[ $LowCov == 1 ]]; then $DIR/bin/05_ObtainSegStats.sh -d $Sdl -D $Sdu -P $Ped -m $mer -L ; else $DIR/bin/05_ObtainSegStats.sh -d $Sdl -D $Sdu -P $Ped -m $mer ; fi && 
if [[ $kin == 1 ]]; then echo -e "\n\e[31mRunning Kmer kinship\e[0m" ; $DIR/bin/05b_KmerKinship.sh -P $Ped -m $mer ; fi &&
if [[ -n $Max ]]; then echo -e "\n\e[31mDownsampling markers\e[0m" ; $DIR/bin/05c_MarkerReduction.sh -P $Ped -m $mer -U $Max ; fi &&
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
