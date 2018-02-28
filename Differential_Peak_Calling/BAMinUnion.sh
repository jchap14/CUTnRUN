#!/bin/bash
##### Trim FASTQs for quality & adaptors. FASTQC output generated also

##### for x in `/bin/ls *.bam` ; do bash BAMinUnion.sh $x ; done
##### bash BAMinUnion.sh $BAM $BED

## add modules
# module add trim_galore

## define variables
BAM=`echo $1`
NAME=`basename $1 .bam`
BED=`echo $2`
# BED="caMEK5_vs_caGFP.K27me3.union.bed"

## write a tempscript to be looped over
cat > $NAME.tempscript.sh << EOF
#!/bin/bash
#$ -N $NAME.BAMtest
#$ -j y
#$ -cwd
#$ -V
#$ -l h_vmem=4G
#$ -pe shm 1
#$ -l h_rt=5:59:00
#$ -l s_rt=5:59:00

##### Commands:
##### test % of BAMs in union peaks for each file 
bedtools intersect -abam $BAM -b $BED > $NAME.union.bam
EOF

## qsub then remove the tempscript
qsub $NAME.tempscript.sh 
sleep 1
rm $NAME.tempscript.sh
