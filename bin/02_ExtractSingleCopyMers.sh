#!/bin/bash -l
#################################################
#	A shell script to generate and plot histograms of parental hashes
#	It will try to calculate peaks, if the user does not define them in the pedigree file, though this may be error prone.
#	Finally, it extracts k-mers it estimates to be single copy
#################################################
Name=$(basename $0)
usage="${Name}; [-h] [-P] [-m] -- A script to obtain single copy k-mers from parental JELLYFISH hashes (AFLAP 2/6).
Options
	-h show this help message
	-P Pedigree file, required. See AFLAP README for more information.
	-m K-mer size. Optional. Default [31]
Temporary files will be output to AFLAP_tmp/02."
#Option block
while getopts ':h:P:m:' option; do
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

DIR=$(dirname $0)
#Check Info from 01_JELLYFISH.sh available.
if [[ -e AFLAP_tmp/01/LA.txt ]]
then
echo -e "\nLinkage analysis to be run using markers derived from:"
cat AFLAP_tmp/01/LA.txt
else
echo -e "\n Could not find output of previous scripts. Please rerun 01_JELLYFISH.sh."
exit 1
fi

#Check JF
JF=$(jellyfish --version)
if [[ $JF =~ ^jellyfish ]]; then echo -e "\n$JF detected" >&2 ; else echo -e "\njellyfish not detected, please modify your PATH" >&2 ; exit 1 ; fi

#Make new tmp for second script.
mkdir -p AFLAP_Results/Plots
mkdir -p AFLAP_tmp/02/ParentalHisto
echo "#Boundaries to follow" > AFLAP_tmp/02/Boundaries.txt
#Obtain Histograms.
echo -e "\nGenerating histograms for F0 to undergo linkage analysis."
F1=$(awk '$2 == 1 {print $4, $5}' AFLAP_tmp/01/Crosses.txt | sort -u | wc -l)
F2=$(awk '$2 == 2 {print $4, $5}' AFLAP_tmp/01/Crosses.txt | sort -u | wc -l)
for g in `cat AFLAP_tmp/01/LA.txt`
  do
	if [[ -e AFLAP_tmp/02/ParentalHisto/$g.${mer}.histo ]] #1
	then
	#Use existing results.
	echo -e "Histogram for $g ${mer}-mer detected. Will use those results.\nNot correct? Cancel and delete AFLAP_tmp/02/ParentalHisto/$g.${mer}.histo, or rerun in clean directory"
	else
	#Generate histogram
	jellyfish histo AFLAP_tmp/01/ParentalCounts/$g.jf${mer} > AFLAP_tmp/02/ParentalHisto/$g.${mer}.histo
	echo -e "Histogram for $g generated"
	fi #1
  #Check to see if Coverage boundaries provided by the user?
  BC=$(awk -v var="$g" '$1 == var && NF == 5 {print $1, $4, $5}' $Ped | head -n 1 | wc -l)
	if [[ $BC == 1 ]] #1
	then
	echo -e "User has supplied lower and upper boundaries for $g ${mer}-mer cut-offs.\nOn to extraction"
	Lo=$(awk -v var="$g" '$1 == var && NF == 5 {print $4}' $Ped | head -n 1)
	Up=$(awk -v var="$g" '$1 == var && NF == 5 {print $5}' $Ped | head -n 1)
	else
	#Super hacky peak finding and boundary setting follows!
	echo -e "User has not supplied boundaries for $g ${mer}-mer cut-offs. Automatic calculation underway\nNote results may not be accurate and plotted figures should be investigated.\n\tIf the F0 are under 20x coverage, reccomend run is cancelled and user defined boundaries be set"
	Peak=$(tail -n+20 AFLAP_tmp/02/ParentalHisto/$g.${mer}.histo | sort -nrk2,2 | head -n 1 | awk '{print $1}')
	echo "Peak found at ${Peak}x. Trying to determing if homozygous or heterozygous."
	HomTest=$(printf %.0f $(echo $Peak*1.5 | bc -l))
	HomPeak=$(tail -n+$HomTest AFLAP_tmp/02/ParentalHisto/$g.${mer}.histo | sort -nrk2,2 | head -n 1 | awk '{print $1}')
	HetTest=$(printf %.0f $(echo $Peak*0.75 | bc -l))
	HetPeak=$(head -n $HetTest AFLAP_tmp/02/ParentalHisto/$g.${mer}.histo | tail -n+20 | sort -nrk2,2 | head -n 1 | awk '{print $1}')
	HomLo=$(printf %.0f $(echo $Peak*1.9 | bc -l))
	HomUp=$(printf %.0f $(echo $Peak*2.1 | bc -l))
	HetLo=$(printf %.0f $(echo $Peak*0.4 | bc -l))
	HetUp=$(printf %.0f $(echo $Peak*0.6 | bc -l))
	Hom=0
	Het=0
		if (( $HomPeak > $HomLo && $HomTest < $HomUp )) #2
		then
		echo "Additional Hom peak detected at $HomPeak"
		Hom=1
		fi #2
		if (( $HetPeak > $HetLo && $HetPeak < $HetUp )) #2
		then
		echo "Additional Het peak detected at $HetPeak"
		Het=1
		fi #2
		if [[ $Hom == 1 && $Het != 1 ]] #2
		then
		echo -e "\nHeterozygous peak =  $Peak\nHomozygous peak = $HomPeak"
			if (( $F1 > 0 )) #3
			then
			echo "Analyzing F1 population so boundaries around the heterozygous peak will be used."
			Lo=$(printf %.0f $(echo $Peak*0.75 | bc -l))
			Up=$(printf %.0f $(echo $Peak*1.5 | bc -l))
			elif (( $F2 > 0 ))
			then
			echo "Analyzing F2 population so boundaries around the homozygous peak will be used."
			Lo=$(printf %.0f $(echo $HomPeak*0.75 | bc -l))
                        Up=$(printf %.0f $(echo $HomPeak*1.5 | bc -l))
			fi #3
		elif [[ $Hom != 1 && $Het == 1 ]]
		then
		echo -e "\nHeterozygous peak = $HetPeak\nHomozygous peak = $Peak"
			if (( $F1 > 0 )) #3
                        then
                        echo "Analyzing F1 population so boundaries around the heterozygous peak will be used."
                        Lo=$(printf %.0f $(echo $HetPeak*0.75 | bc -l))
                        Up=$(printf %.0f $(echo $HetPeak*1.5 | bc -l))
                        elif (( $F2 > 0 ))
                        then
                        echo "Analyzing F2 population so boundaries around the homozygous peak will be used."
                        Lo=$(printf %.0f $(echo $Peak*0.75 | bc -l))
                        Up=$(printf %.0f $(echo $Peak*1.5 | bc -l))
                        fi #3
		elif [[ $Hom == 1 && $Het == 1 ]]
		then
		echo -e "\nLooks like three peaks have been found. That's going to make things tricky!\nPeak 1 (Het) = $HetPeak\nPeak 2 (Het) = $Peak\nPeak 3 (Hom) = $HomPeak"
			if (( $F1 > 0 )) #3
                        then
                        echo "Analyzing F1 population so boundaries around the two heterozygous peaks will be used."
                        Lo=$(printf %.0f $(echo $HetPeak*0.75 | bc -l))
                        Up=$(printf %.0f $(echo $Peak*1.5 | bc -l))
                        elif (( $F2 > 0 ))
                        then
                        echo "Analyzing F2 population so boundaries around the homozygous peak will be used, however it seems like there is a large amount of heterozygosity"
                        Lo=$(printf %.0f $(echo $HomPeak*0.75 | bc -l))
                        Up=$(printf %.0f $(echo $HomPeak*1.5 | bc -l))
                        fi #3
		elif [[ $Hom != 1 && $Het != 1 ]]
		then
		echo -e "\nOnly one peak found. This will be used to obtain markers."
		Lo=$(printf %.0f $(echo $Peak*0.75 | bc -l))
		Up=$(printf %.0f $(echo $Peak*1.5 | bc -l))
 		fi #2
	echo -e "Peak boundary estimation complete. They are visible on AFLAP_Results/Plots/${g}_m${mer}_L${Lo}_U${Up}_histo.png\n"
	fi #1
  echo "Lower boundary for $g set to $Lo, upper boundary to $Up"
  echo -e "$g\t$Lo\t$Up" >> AFLAP_tmp/02/Boundaries.txt
  if [[ -e AFLAP_tmp/02/ParentalHisto/${g}_m${mer}_L${Lo}_U${Up}.fa ]]
  then 
  echo -e "\n$g ${mer}-mer previously extracted between $Lo and $Up. Delete ${g}_m${mer}_L${Lo}_U${Up}.fa to rebuild." 
  else
  jellyfish dump -U $Up -L $Lo -o AFLAP_tmp/02/ParentalHisto/${g}_m${mer}_L${Lo}_U${Up}.fa AFLAP_tmp/01/ParentalCounts/$g.jf${mer}
  Kco=$(grep -c '^>'  AFLAP_tmp/02/ParentalHisto/${g}_m${mer}_L${Lo}_U${Up}.fa) 
  echo "$Kco ${mer}-mers extracted from $g"
  Rscript $DIR/HistoPlot.R AFLAP_tmp/02/ParentalHisto/${g}.${mer}.histo $Lo $Up AFLAP_Results/Plots/${g}_m${mer}_L${Lo}_U${Up}_histo.png
  fi
  done
exit
