#!/bin/bash -l
#Usage Statement

Name=$(basename $0)
usage="${Name}; [-h] [-P] [-t] [-m] -- A script to obtain JELLYFISH hashes to be used when running AFLAP.

Options
	-h show this help message
	-P Pedigree file, required. See AFLAP README for more information.
	-m K-mer size. Optional. Default [31]
	-t Threads for JELLYFISH counting. Optional. Default [4]

Temporary files will be output to AFLAP_tmp/01."

while getopts ':hP:t:m:' option; do
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
#Set defaults if not specified.
if [[ -z $thread ]]; then
echo "Threads not specified, will proceed with default [4]"
thread=4
fi

if [[ -z $mer ]]; then
echo "mer size not specified, will proceed with default [31]"
mer=31
fi

#0. Dependency Check
#JELLYFISH
JF=$(jellyfish --version)
if [[ $JF =~ ^jellyfish ]]; then echo "$JF detected"; else echo "jellyfish not detected, please modify your PATH" ; exit 1 ; fi

#Make a tmp directory
mkdir -p AFLAP_tmp/01
mkdir -p AFLAP_Intermediate

#Pedigree file check
#Check labels
Flab=$(awk '$2 > 2' $Ped | wc -l)
if [[ $Flab > 0 ]]; then echo "Pedigree file specifies individuals which are not F0, F1 or F2. AFLAP has not been validated on this type of data so will terminate." ; exit 1 ; fi
Ilab=$(awk '{print $1, $2}' $Ped | sort -u | awk '{print $1}' | sort | uniq -c | awk '$1 > 1' | wc -l)
if [[ $Ilab > 0 ]]; then echo "Pedigree file specifies the same individual label over multiple generations. This is incompatable with AFLAP so will terminate." 
awk '{print $1, $2}' $Ped | sort | awk '{print $1}' | sort | uniq -c | awk '$1 > 1'
exit 1
fi

#Calculate number of parents in pedigree file:
awk '$2 == 0 {print $1}' $Ped | sort -u > AFLAP_tmp/01/F0.txt
Pcou=$(wc -l AFLAP_tmp/01/F0.txt | awk '{print $1}')

if [[ $Pcou == 2 ]]
then
echo -e "$Pcou parents detected. Simple! Beginning k-mer counting"
sort -uk1,1 $Ped | awk '$2 > 0 {print $2, $4, $5}' | sort | uniq -c > AFLAP_tmp/01/Crosses.txt
elif [[ $Pcou > 2 ]]
then
echo -e "$Pcou parents detected. Not so simple! Determining crosses"
sort -uk1,1 $Ped | awk '$2 > 0 {print $2, $4, $5}' | sort | uniq -c > AFLAP_tmp/01/Crosses.txt
CrossCou=$(cat AFLAP_tmp/01/Crosses.txt | wc -l)
	if [[ $Crou > 4 ]]
	then
	echo "$CrossCou crosses identified in pedigree file. This may take a long time to process. Crosses include:"
	sort -nrk AFLAP_tmp/01/Crosses.txt | head -n 5 | awk '{print $3, $4}'
	else
	echo "$CrossCou crosses identified:"
	cat AFLAP_tmp/01/Crosses.txt | awk '{print $3, $4}'
	fi
elif [[ $Pcou < 2 ]]
then
echo -e "$Pcou parents detected. Two or more F0 should be specified in the pedigree file to run AFLAP. Exiting."
exit 1
fi

#Check pedigree stucture.
#Currently only tested on F1 and F2 independently. Has not been tested on both simultaneously.
F1=$(awk '$2 == 1 {print $3, $4}' AFLAP_tmp/01/Crosses.txt | sort -u | wc -l)
F2=$(awk '$2 == 2 {print $3, $4}' AFLAP_tmp/01/Crosses.txt | sort -u | wc -l)
if [[ $F1 > 0 && $F2 > 0 ]]
then
echo "both F1 and F2 populations detected. AFLAP has not been validated on this type of data so will terminate. Please try to rerun specifying either type of population."
exit 1
elif [[ $F1 > 0 && $F2 == 0 ]]
then
echo "$F1 F1 cross(es) identified."
awk '$2 == 1 {print $3 "x" $4}' AFLAP_tmp/01/Crosses.txt > AFLAP_tmp/01/ParentsToCompare.txt
head AFLAP_tmp/01/ParentsToCompare.txt
elif [[ $F1 == 0 && $F2 > 0 ]]
then
echo "$F2 F2 cross(es) identified."
awk '$2 == 1 {print $3 "x" $4}' AFLAP_tmp/01/Crosses.txt > AFLAP_tmp/01/ParentsToCompare.txt
head AFLAP_tmp/01/ParentsToCompare.txt
fi

#Identify which parents the map is to be constructed of.
#Mostly this information will be used in later scripts.
awk '$2 == 0 && $4 ~ /NA/ {print $1}' $Ped | sort -u > AFLAP_tmp/01/NoLA.txt
awk '$2 == 0 {print $1}' $Ped | sort -u | join -v 1 - AFLAP_tmp/01/NoLA.txt > AFLAP_tmp/01/LA.txt
LMB=$(cat AFLAP_tmp/01/LA.txt | wc -l)
NLM=$(cat AFLAP_tmp/01/NoLA.txt | wc -l)
if [[ $LMB == 0 ]]
then
echo -e "\nHmm it seems the user has requested no linkage maps be built with parental derived markers.\nI reccomend you check the pedgree file, atleast one F0 should have no NA in fields 4 and 5. \nThey can be blank for automatic coverage calculation, or contain coverage cut-off boundaries."
exit 1
else
echo -e "\nWe are in buisness! Linkage analysis will be performed on markers derived from $LMB F0:"
cat AFLAP_tmp/01/LA.txt
        if [[ $NLM >1 ]]
        then
        echo -e "\nLinkage analysis will not be performed on these F0."
	cat AFLAP_tmp/01/NoLA.txt
        echo -e "\nDoesn't sound right? Make sure 'NA' does not occur in fields 4 and 5 for these F0 in $Ped"
	fi
echo -e "Moving onto k-mer counting\n"
fi
#K-mer counting.
#Parents.
mkdir -p AFLAP_Intermediate/ParentalCounts
PLA=$(cat AFLAP_tmp/01/F0.txt | wc -l)
echo -e "Begining k-mer counting for $PLA parents"
for g in `cat AFLAP_tmp/01/F0.txt`
  do
#Check if has has already been created
	if [[ -e AFLAP_Intermediate/ParentalCounts/$g.jf${mer} ]]
	then
#If so, skip.
	echo -e "$mer hash detected for $g. Skipping. \n\tNot correct? Cancel and delete AFLAP_Intermediate/ParentalCounts/$g.fj${mer}, or run from clean directory"
	else
#If not, generate it.
	echo -e "Begining k-mer counting for $g"
	Reads=$(awk -v var="$g" '$1 == var {print $3}' $Ped | tr '\n' ' ')
	jellyfish count -s 1G -t $thread -m $mer -C -o AFLAP_Intermediate/ParentalCounts/$g.jf${mer} <(zcat $Reads) &&
	fi
  done
echo -e "\nParental K-mer counting complete!\nOn to the progeny!\n"
#Progeny.
#Would be good to parallelize.
mkdir -p AFLAP_Intermediate/ProgCounts
awk '$2 != 0 {print $1}' $Ped | sort -u > AFLAP_tmp/01/Prog.txt
ProgCou=$(cat AFLAP_tmp/01/Prog.txt | wc -l)
echo -e "Begining k-mer counting for $ProgCou progeny"
for g in `cat AFLAP_tmp/01/Prog.txt`
  do
#Check if has has already been created
        if [[ -e AFLAP_Intermediate/ProgCounts/$g.jf${mer} ]]
        then
#If so, skip.
        echo -e "$mer hash detected for $g. Skipping. \n\tNot correct? Cancel and delete AFLAP_Intermediate/ProgCounts/$g.fj${mer}, or run from clean directory"
        else
#If not, generate it.
	echo -e "Begining k-mer counting for $g"
	Reads=$(awk -v var="$g" '$1 == var {print $3}' $Ped | tr '\n' ' ')
	jellyfish count -s 1G -t $thread -m $mer -C -o AFLAP_Intermediate/ProgCounts/$g.jf${mer} <(zcat $Reads) &&
	fi
  done
echo -e "\nk-mer counting done!"
exit
#for g in `seq 1 1 $Pcou` ; do echo $g ; done
#Were user specified boundaries provided:


