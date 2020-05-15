#!/bin/bash -l
#Usage Statement

while getopts ':hP:' option; do
        case "$option" in
                h)  echo "$usage"
                         exit
                         ;;
                P)  Ped=$OPTARG
                         ;;
                \?) printf "illegal option: -%s\n\n" "$OPTARG" >&2
#                    echo "$usage"
                    exit 1
                         ;;
        esac
done

echo "Beginning AFLAP"

#0. Dependency Check
#JELLYFISH
JF=$(jellyfish --version)
if [[ $JF =~ ^jellyfish ]]; then echo "$JF detected"; else echo "jellyfish not detected, please modify your PATH" ; exit 1 ; fi

#Make a tmp directory
mkdir -p 00_AFLAPtmp

#Pedigree file check
#Check labels
Flab=$(awk '$2 > 2' $Ped | wc -l)
if [[ $Flab > 0 ]]; then echo "Pedigree file specifies individuals which are not F0, F1 ir F2. AFLAP has not been validated on this type of data so will terminate." ; rm -r 00_AFLAPtmp ; exit 1 ; fi
Ilab=$(awk '{print $1, $2' $Ped | sort | awk '{print $1}' | sort | uniq -c | awk '$2 > 1' | wc -l)
if [[ $Ilab > 0 ]]; then echo "Pedigree file specifies the same individual label over multiple generations. This is incompatable with AFLAP so will terminate." 
awk '{print $1, $2' $Ped | sort | awk '{print $1}' | sort | uniq -c | awk '$2 > 1' 
rm -r 00_AFLAPtmp 
exit 1 
fi


#Calculate number of parents in pedigree file:
awk '$2 == 0 {print $1}' $Ped | sort -u > 00_AFLAPtmp/F0.txt
Pcou=$(wc -l 00_AFLAPtmp/F0.txt | awk '{print $1}')

if [[ $Pcou == 2 ]]
then
echo -e "$Pcou parents detected. Simple! Beginning k-mer counting"
else if [[ $Pcou > 2 ]]
then
echo -e "$Pcou parents detected. Not so simple! Determining crosses"

#Split into seperate loop from here
#Should add condition that there is at lease one common parent

sort -uk1,1 $Ped | awk '$2 > 0 {print $2, $4, $5}' | sort | uniq -c > 00_AFLAPtmp/Crosses.txt
F1=$(awk '$2 == 1 {print $3, $4}' 00_AFLAPtmp/Crosses.txt | sort -u | wc -l)
F2=$(awk '$2 == 2 {print $3, $4}' 00_AFLAPtmp/Crosses.txt | sort -u | wc -l)
if [[ F1 > 0 && F2 > 0 ]]
then
echo "both F1 and F2 populations detected. AFLAP has not been validated on this type of data so will terminate. Please try to rerun specifying either type of population."
rm -r 00_AFLAPtmp
exit 1
elif [[ F1 > 0 && F2 == 0]]
then
echo "$F1 F1 crosses identified. Will use this information when obtaining parental specific markers."
awk '$2 == 1 {print $3, $4}' tmp/Crosses.txt | sort -u | uniq -c > 00_AFLAPtmp/ParentsToCompare.txt
head 00_AFLAPtmp/ParentsToCompare.txt
elif [[ F1 == 0 && F2 > 0]]
then
echo "$F2 F2 crosses identified. Will use this information when obtaining parental specific markers."
awk '$2 == 1 {print $3, $4}' tmp/Crosses.txt | sort -u | uniq -c > 00_AFLAPtmp/ParentsToCompare.txt
head 00_AFLAPtmp/ParentsToCompare.txt
fi

