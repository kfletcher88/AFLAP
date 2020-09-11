#!/bin/bash -l

while getopts ':hP:m:L:T:' option; do
        case "$option" in
                h)  echo "$usage"
                         exit
                         ;;
                P)  Ped=$OPTARG
                         ;;
                m)  mer=$OPTARG
                         ;;
		L)  LOD=$OPTARG
			;;
		T)  Threads=$OPTARG
			;;
                \?) printf "illegal option: -%s\n\n" "$OPTARG" >&2
                    echo "$usage"
                    exit 1
                         ;;
        esac
done

#Check if options aren't set.

if [[ -e AFLAP_tmp/01/LA.txt && -e AFLAP_tmp/02/Boundaries.txt ]]
then
echo -e "\nIntermediate files detected"
else
echo -e "\n Could not find output of previous scripts. Please rerun full pipeline."
exit 1
fi

DIR=$(dirname $0)
mkdir -p AFLAP_Results/LOD$LOD

if [[ ! -e $DIR/../ThirdParty/LepMap3/bin/SeparateChromosomes2.class ]]
then
echo -e "\nLepMap3 modules not detected in $DIR/../ThirdParty/LepMap3/bin/\n\nPlease run $DIR/LepMap3Installer.sh to install LepMap3 here.\n\tAlternatively, run LepMap3 by providing the installation directory with option -d to perl scripts: \n\t\t$DIR/LepMap3SeperateChromosomes.pl and \n\t\t$DIR/LepMap3OrderMarkers.pl\nExiting"
exit 1
fi

#Set variables to obtain correct genotype table.
for g in `cat AFLAP_tmp/01/LA.txt`
do
        if [[ -e AFLAP_tmp/03/${g}_CrossedTo.txt ]]
                then
                P0=$(cat AFLAP_tmp/03/${g}_CrossedTo.txt | tr '\n' '_' | sed 's/_$//')
                Lo=$(awk -v var="$g" '$1 == var {print $2}' AFLAP_tmp/02/Boundaries.txt)
                Up=$(awk -v var="$g" '$1 == var {print $3}' AFLAP_tmp/02/Boundaries.txt)
                else
                echo -e "Intermediate file is missing, please rerun AFLAP (s0603P1CT)"
                exit 1
        fi
	if [[ ! -f AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}.ForLepMap3.tsv ]]; then echo "Genotype table is missing. Please rerun AFLAP.sh to generate" ; exit 1 
	else
	echo -e "Initiating LepMap3 for $g\n"
	if [[ -f AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.txt ]]; then echo "Previous results detected. Skipping Seperate Chromosomes"
	else
	java -cp $DIR/../ThirdParty/LepMap3/bin/ SeparateChromosomes2 data=AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}.ForLepMap3.tsv lodLimit=$LOD numThreads=$Threads > AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.txt 2>AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.stderr
	awk 'NR > 1' AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.txt | sort -n | uniq -c > AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.fre
	fi
	Mcou=$(awk '{sum += $1}END{print sum}' AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.fre)
	LGcou=$(awk -v Mcou=$Mcou '$2 != 0 && $1/Mcou >= 0.01' AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.fre | wc -l)
	echo -e "$LGcou Linkage groups detected, containing a minimum of 1% of the markers"
		if (( $Threads > $LGcou ))
		then echo -e "Running linkage group ordering in parallel as number of threads exceeds number of linkage groups (Maximum efficiency!)"
		for L in `awk -v Mcou=$Mcou '$2 != 0 && $1/Mcou >= 0.01 {print $2}' AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.fre` ; do if [[ ! -f AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.LG$L.txt ]]; then java -cp $DIR/../ThirdParty/LepMap3/bin/ OrderMarkers2 useMorgan=1 numMergeIterations=20 chromosome=$L data=AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}.ForLepMap3.tsv map=AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.txt > AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.LG$L.txt 2> AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.LG$L.stderr & fi ; done
		wait
		else
		echo -e "Running linkage group ordering in series as number of threads does not exceed number of linkage groups"
		for L in `awk -v Mcou=$Mcou '$2 != 0 && $1/Mcou >= 0.01 {print $2}' AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.fre` ; do if [[ ! -f AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.LG$L.txt ]]; then java -cp $DIR/../ThirdParty/LepMap3/bin/ OrderMarkers2 useMorgan=1 numMergeIterations=20 chromosome=$L data=AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}.ForLepMap3.tsv map=AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.txt > AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.LG$L.txt 2> AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.LG$L.stderr ; fi ; done
		fi
		echo "Linkage group ordering complete"
	fi
#Use Pedigree table to see if M/F map built. Extract information. Build marker fasta.
	Z=$(awk -v var=$g '$4 == var || $5 == var {if ($4 == var) print 4; else if ($5 == var) print 5}' $Ped | head -n 1)
	if (( Z == 4 )); then
		cat AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.LG*.txt | awk -v OFS='\t' '$2 ~ /LG/ {lg = $4; next} $1 !~ /\#/ {print lg, $2, $1}' | sort -k3,3 | \
		join -1 3 - <(cut -f1,2 AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}.ForLepMap3.tsv | awk '(NR>6){print ++i, $0}' | sort -k1,1) | awk -v OFS='\t' '{print $4, $2, $3, $5}' > AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD$LOD.txt
	elif (( Z == 5 )); then
		cat AFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.LG*.txt | awk -v OFS='\t' '$2 ~ /LG/ {lg = $4; next} $1 !~ /\#/ {print lg, $3, $1}' | sort -k3,3 | \
		join -1 3 - <(cut -f1,2 AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}.ForLepMap3.tsv | awk '(NR>6){print ++i, $0}' | sort -k1,1) | awk -v OFS='\t' '{print $4, $2, $3, $5}' > AFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD$LOD.txt
	fi
echo -e "Genetic map for $g:\n\tAFLAP_Results/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD$LOD.txt\nLepMap3 output:\n\tAFLAP_Results/LOD${LOD}\nLepMap3 Map file:\n\tAFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}.txt\n\n"
echo -e "Appending to log file"
if [[ ! -f AFLAP_RUN.LOG ]]; then echo -e "Date\tMappedParent\tOtherParents\tmerLen\tUpperBoundary\tLowerBoundary\tMinLOD\tLepMap3Outputs\tLinkageGroupCount" >> AFLAP_RUN.LOG ; fi
paste <(date) <(echo -e "$g\t$P0\t$mer\t$Lo\t$Up\t$LOD\tAFLAP_Results/LOD${LOD}/${g}_m${mer}_L${Lo}_U${Up}_${P0}.LOD${LOD}*\t$LGcou") >> AFLAP_RUN.LOG
done
echo -e "AFLAP Finished!!!"
exit
