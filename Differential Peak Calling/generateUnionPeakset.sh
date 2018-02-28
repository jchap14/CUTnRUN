#!/bin/bash

## run command: bash generateUnionPeakset.sh $.metadata.csv

##### concatenate all peaksets
METAFILE=`echo $1` #'caMEK5_vs_caGFP.K27me3.metadata.csv'
NAME=`basename $METAFILE .metadata.csv`
mac2unix $METAFILE
COLNUM=`head -1 $METAFILE| tr ',' '\n' | nl | grep 'Peaks' | cut -f 1`
cat $METAFILE | cut -d "," -f $COLNUM | tail -n+2 > PEAKSaddress.txt

##### optional: add a window to each side (e.g. 100bp)

##### union set by bedSort + bedMerge
cat `cat PEAKSaddress.txt | tr '\n' ' '` | cut -f1,2,3 > tempA.bed
sortBed -i tempA.bed | mergeBed -i stdin > $NAME.union.bed
rm tempA.bed
