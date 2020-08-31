#!/bin/bash -l
#Make directory, surpress error.
mkdir -p ReadsEBI
#Change into the directory for download
cd ReadsEBI
#Set variabled to access EBI ftp
Sample=$1
S1=${Sample: 0:6}
S2=${Sample: -1}
#Download
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${S1}/00${S2}/${Sample}/${Sample}_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${S1}/00${S2}/${Sample}/${Sample}_2.fastq.gz

