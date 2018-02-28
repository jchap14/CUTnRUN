#!/bin/bash
##### Use this script to remove GFP control peaks from CUT&RUN peaks

##### submission command
## bash GFPFilterPeaks.sh [$NAME.filt.Peaks.bed] [CONTROLFILE.]

##### load required modules (comment for local use)
module add bedtools

##### set variables
NAME=`basename $1 .filt.Peaks.bed`
PEAKfile=`echo $1`
CONTROLFILE=`echo $2`

##### remove peaks that intersect with the Control GFP file
## Only report those entries in A that have no overlap in B
bedtools intersect -v -a $PEAKfile -b $CONTROLFILE > $NAME.enrich.Peaks.bed
