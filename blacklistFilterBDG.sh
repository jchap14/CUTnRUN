#!/bin/bash
##### Use this script to remove highly mappable blacklisted regions from CUT&RUN peaks

##### submission command
## for x in `/bin/ls *.spikenorm.bdg` ; do bash blacklistFilterBDG.sh $x; done

##### load required modules
module add bedtools

##### set variables
NAME=`basename $1 .spikenorm.bdg`
BDGfile=`echo $1`
BLACKLIST=/srv/gsfs0/projects/snyder/chappell/Annotations/UCSC-hg19/hg19EncodeMapabilityBlacklist.bed

##### create a tempscript for queue sub
cat > $NAME.tempscript.sh << EOF
#!/bin/bash
#$ -N $NAME.filtBDG
#$ -j y
#$ -V
#$ -cwd
#$ -l h_vmem=5G
#$ -pe shm 1

##### remove peaks that intersect with the blacklist.bed
## Only report those entries in A that have no overlap in B
bedtools intersect -v -a $BDGfile -b $BLACKLIST > $NAME.filt.bdg
EOF

## qsub then remove the tempscript
qsub $NAME.tempscript.sh 
sleep 1
rm $NAME.tempscript.sh
