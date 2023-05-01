#!/bin/bash -l
for g in `cat Accession.txt` ; do ../bin/00_ReadDload.sh $g ; done 
