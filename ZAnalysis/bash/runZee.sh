#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Illegal number of parameters"
    echo "<inputFileName> <Tag> <isMC>"
    exit
fi

inDir=$1
tag=$2
isMC=$3

for i in $(seq 0 14);
do
    ./bin/Z_EE_Channel.exe $inDir $tag $isMC $i 15 > unmergedOutputs/job_$Tag_$isMC_$i.log 2>&1 & 
done
