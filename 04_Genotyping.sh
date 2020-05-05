#!/bin/bash -l

#Usage Statement

while getopts ':hd:1:2:' option; do
        case "$option" in
                h)  echo "$usage"
                         exit
                         ;;
                d)  dir=$OPTARG
                         ;;
                1)  M61=$OPTARG
                         ;;
                2)  M62=$OPTARG
                         ;;
                \?) printf "illegal option: -%s\n\n" "$OPTARG" >&2
                    echo "$usage"
                    exit 1
                         ;;
        esac
done

#0. Dependency Check
#a JELLYFISH
JF=$(jellyfish --version)
if [[ $JF =~ ^jellyfish ]]; then echo "$JF detected"; else echo "jellyfish not detected, please modify your PATH" ; exit 1 ; fi

#P labels could be taken from a Pedigree file
mkdir 04_AFLAP_P1_61bp
mkdir 04_AFLAP_P1_62bp
mkdir 04_AFLAP_P2_61bp
mkdir 04_AFLAP_P2_62bp

#Loop
for g in $dir ; do o=$(basename $g | sed 's/$/_61bp.count/') ; jellyfish query -s $M61 $g > 04_AFLAP_P1_61bp/$o ; done
for g in $dir ; do o=$(basename $g | sed 's/$/_62bp.count/') ; jellyfish query -s $M62 $g > 04_AFLAP_P1_62bp/$o ; done
