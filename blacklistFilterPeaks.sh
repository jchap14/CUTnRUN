#!/bin/bash
##### Use this script to remove highly mappable blacklisted regions from CUT&RUN peaks

##### submission command
## for x in `/bin/ls *.narrowPeak` ; do bash blacklistFilterPeaks.sh $x; done
## for x in `/bin/ls *.peaks.bed` ; do bash blacklistFilterPeaks.sh $x; done

##### load required modules (uncomment for SCG use)
# module add bedtools

##### set variables
# NAME=`basename $1 .narrowPeak`
NAME=`basename $1 .peaks.bed`
PEAKfile=`echo $1`
## uncomment for SCG use 
# BLACKLIST=/srv/gsfs0/projects/snyder/chappell/Annotations/UCSC-hg19/hg19EncodeMapabilityBlacklist.bed
## uncomment for local use 
BLACKLIST="/Users/jchap12/Google\ Drive/BIOINFORMATICS/Annotations/hg19/hg19EncodeMapabilityBlacklist.bed"

##### create a tempscript for queue sub
cat > $NAME.tempscript.sh << EOF
#!/bin/bash
#$ -N $NAME.filtBlacklist
#$ -j y
#$ -V
#$ -cwd
#$ -l h_vmem=5G
#$ -pe shm 1

##### remove peaks that intersect with the blacklist.bed
## Only report those entries in A that have no overlap in B

bedtools intersect -v -a $PEAKfile -b $BLACKLIST > $NAME.filt.Peaks.bed
EOF

## qsub on scg
# qsub $NAME.tempscript.sh
## bash locally
bash $NAME.tempscript.sh

sleep 1
## remove the tempscript
rm $NAME.tempscript.sh
