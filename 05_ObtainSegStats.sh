#!/bin/bash -l
#Obtain Segregation Histograms
awk '{ for(i=2; i<=NF;i++) j+=$i; print j; j=0 }' AFLAP_61bp_GT.tsv | sort -n | uniq -c > AFLAP_61bp_GT.hist
awk '{ for(i=2; i<=NF;i++) j+=$i; print j; j=0 }' AFLAP_62bp_GT.tsv | sort -n | uniq -c > AFLAP_62bp_GT.hist
cat AFLAP_61bp_GT.tsv AFLAP_62bp_GT.tsv | awk '{ for(i=2; i<=NF;i++) j+=$i; print j; j=0 }' | sort -n | uniq -c > AFLAP_AllMarkers_GT.hist

