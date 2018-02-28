#!/bin/bash
##### generate normalized yeast coverage bedgraph for CUT&RUN

##### for x in `/bin/ls *.yeast.bedpe` ; do bash CUT_RUN_spike_norm_coverage_bedgraph.sh $x; done

## load required modules
module add bedtools/2.21.0

## define variables
NAME=`basename $1 .yeast.bedpe`
GENOME_FILE="/srv/gsfs0/projects/snyder/chappell/Annotations/UCSC-hg19/hg19.chrom.sizes"

## calculating yeast MAPPED_NUM
echo "calculating yeast MAPPED_NUM:"
MAPPED_NUM=`wc -l $NAME.yeast.bedpe | cut -f1 -d ' '`
echo $MAPPED_NUM

## calculating FACTOR
echo "calculating FACTOR"
FACTOR=`awk "BEGIN {print 10000/\$MAPPED_NUM}"`
echo $FACTOR

##### write a tempscript to be looped over
cat > $NAME.tempscript.sh << EOF
#!/bin/bash
#$ -N $NAME.spikeNorm
#$ -j y
#$ -cwd
#$ -V
#$ -l h_vmem=20G
#$ -pe shm 1
#$ -l h_rt=0:59:00
#$ -l s_rt=0:59:00

## run commands
echo "sorting bedpe"
bedtools sort -i $NAME.bedpe > $NAME.sorted.bedpe
echo "calculating coveraage from bedpe"
bedtools genomecov -bga -scale $FACTOR -i $NAME.sorted.bedpe -g $GENOME_FILE > $NAME.spikenorm.bdg
EOF

## qsub then remove the tempscript
qsub $NAME.tempscript.sh 
sleep 1
rm $NAME.tempscript.sh
