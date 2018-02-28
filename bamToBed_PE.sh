#!/bin/bash
##### convert BAMs to bedPE (paired end BAMs)

# for x in `/bin/ls *.bam` ; do bash bamToBed_PE.sh $x; done

## add modules
module add samtools
module add bedtools

## define variables
NAME=`basename $1 .bam`

cat > $NAME.tempscript.sh << EOF
#!/bin/bash
#$ -N $NAME.bamToBed_PE
#$ -j y
#$ -cwd
#$ -V
#$ -l h_vmem=4G
#$ -pe shm 12
#$ -l h_rt=11:59:00
#$ -l s_rt=11:59:00

## run commands
# Sort by read name
samtools sort -n $NAME.bam -o $NAME.sorted.bam
# convert all reads, or...
bedtools bamtobed -i $NAME.sorted.bam -bedpe > $NAME.bedpe
EOF

## qsub then remove the tempscript
qsub $NAME.tempscript.sh
sleep 1
rm $NAME.tempscript.sh
