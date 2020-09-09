#!/bin/bash -l
#Shell script to install LepMap3
#Tested Spetember 9th 2020
#Please visit https://sourceforge.net/projects/lep-map3/
DIR=$(dirname $0)
cd $DIR/../ThirdParty/LepMap3 &&
wget https://sourceforge.net/projects/lep-map3/files/binary%2Bcode.zip &&
unzip binary\+code.zip &&
rm binary\+code.zip &&
mv README LepMap3_README &&
exit
