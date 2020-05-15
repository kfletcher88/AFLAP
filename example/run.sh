#!/bin/bash -l
for g in `cat Accession.txt` ; do ../00_NCBIdload.sh $g ; done 
