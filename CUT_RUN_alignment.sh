##### Align CUT&RUN reads to UCSC hg19 & FlyBase r6.06

##### submit for all FQs in cwd
# for x in `/bin/ls *.trim.R1.fq.gz` ; do bash CUT_RUN_alignment.sh $x; done

## add modules
module add bowtie/2.2.6
# module add samtools ## samtools is in conda env

## define variables
name=`basename $1 .trim.R1.fq.gz`

## write a tempscript to be looped over
cat > $name.tempscript.sh << EOF
#!/bin/bash
#$ -N $name.CR_align
#$ -j y
#$ -cwd
#$ -V
#$ -l h_vmem=6G
#$ -pe shm 4
#$ -l h_rt=11:59:00
#$ -l s_rt=11:59:00

## run commands
##usage: bamToBed [OPTIONS] -i <BAM>
## run bowtie against human genome
echo "Starting Bowtie Alignment against Human"
bowtie2 --local --very-sensitive-local --no-unal --no-mixed \
--no-discordant -I 10 -X 700 --threads 4 \
-x /ifs/scratch/leepc12/pipeline_genome_data/hg19/bowtie2_index/male.hg19.fa \
-1 $1 -2 $name.trim.R2.fq.gz -S $name.sam
# sorts by coordinate -> works for flagstat and idxstats
samtools view -Su $name.sam | samtools sort -o $name.bam
rm $name.sam

## run bowtie against fly genome
echo "Starting Bowtie Alignment against Yeast"
bowtie2 --local --very-sensitive-local --no-unal --no-mixed \
--no-discordant -I 10 -X 700 --threads 4 \
-x /srv/gsfs0/projects/snyder/chappell/Annotations/yeast/S288C \
-1 $1 -2 $name.trim.R2.fq.gz -S $name.yst.sam
samtools view -Su $name.yst.sam | samtools sort -o $name.yst.bam
rm $name.yst.sam
EOF

## qsub then remove the tempscript
qsub $name.tempscript.sh 
sleep 1
rm $name.tempscript.sh
