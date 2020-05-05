#!/bin/bash -l

while getopts ':hj:p:L:U:' option; do
        case "$option" in
                h)  echo "$usage"
                         exit
                         ;;
                j)  Hash=$OPTARG
                         ;;
                U)  Up=$OPTARG
                         ;;
                L)  Lo=$OPTARG
                         ;;
                p)  Out=$OPTARG
                         ;;
                \?) printf "illegal option: -%s\n\n" "$OPTARG" >&2
                    echo "$usage"
                    exit 1
                         ;;
        esac
done


#Check JF
JF=$(jellyfish --version)
if [[ $JF =~ ^jellyfish ]]; then echo "$JF detected" >&2 ; else echo "jellyfish not detected, please modify your PATH" >&2 ; exit 1 ; fi

jellyfish dump -U $Up -L $Lo -o ${Out}_L${Lo}_U${Up}.fa $Hash
 
Kco=$(grep -c '^>' ${Out}_L${Lo}_U${Up}.fa)
echo "$Kco k-mers extracted from $Hash using $JF" >&2 
exit
