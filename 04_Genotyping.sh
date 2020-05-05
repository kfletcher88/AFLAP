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

#Covert to 1s and zeros
#Need a low cov and high cov option
mkdir 04_AFLAP_P1_61bp_GT
mkdir 04_AFLAP_P1_62bp_GT
for g in 04_AFLAP_P1_61bp/* ; do o=$(basename $g) ; awk '{if($2 > 0) print 1; else print 0}' $g > 04_AFLAP_P1_61bp_GT/$o ; done
for g in 04_AFLAP_P1_62bp/* ; do o=$(basename $g) ; awk '{if($2 > 0) print 1; else print 0}' $g > 04_AFLAP_P1_62bp_GT/$o ; done

#Build Genotype table
Mseq=$(ls 04_AFLAP_P1_61bp/* | head -n 1)
paste <(awk '{print $1}' $Mseq) 04_AFLAP_P1_61bp_GT/* > AFLAP_61bp_GT.tsv
Mseq=$(ls 04_AFLAP_P1_62bp/* | head -n 1)
paste <(awk '{print $1}' $Mseq) 04_AFLAP_P1_62bp_GT/* > AFLAP_62bp_GT.tsv

#Obtain Segregation Histograms
awk '{ for(i=2; i<=NF;i++) j+=$i; print j; j=0 }' AFLAP_61bp_GT.tsv | sort -n | uniq -c > AFLAP_61bp_GT.hist
awk '{ for(i=2; i<=NF;i++) j+=$i; print j; j=0 }' AFLAP_62bp_GT.tsv | sort -n | uniq -c > AFLAP_62bp_GT.hist
cat AFLAP_61bp_GT.tsv AFLAP_62bp_GT.tsv | awk '{ for(i=2; i<=NF;i++) j+=$i; print j; j=0 }' | sort -n | uniq -c > AFLAP_AllMarkers_GT.hist
exit
