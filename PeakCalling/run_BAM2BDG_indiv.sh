#!/bin/bash
# This is the submission script for "BAM2BDG_indiv.py"

# submit for a specific bamfile
# bash run_BAM2BDG_indiv.sh $SAMPLE.nmSort.bam

##### submit for all sorted bedgraphs in CWD
# for x in `/bin/ls *.nmSort.bam` ; do bash run_BAM2BDG_indiv.sh $x; done

##### add modules
source activate CUT_n_RUN
module add ucsc_tools/2.7.2


##### specify variables to pass to BAM2BDG_indiv.py
BAMFILE=$1
NAME=`basename $BAMFILE .nmSort.bam`
SPIKEFILE=`echo $NAME.yst.bam`
MERGE_CLOSE_PEAKS="True"
LENGTHS_ANALYSIS=True
LENGTHS_IMAGE=True
SIZE_SELECT_1=True
SIZE_SELECT_2=True
BEDGRAPH=True
CHROM_SIZES="hg19"
BIG_WIG=True
SIZE_MIN_1=20
SIZE_MAX_1=120
SIZE_MIN_2=150
SIZE_MAX_2=710
MULTIPLYING_FACTOR=10000

## create tempscript
cat > $NAME.tempscript.sh << EOF
#!/bin/bash
#$ -N $NAME.BAM2BDG
#$ -S /bin/bash
#$ -j y
#$ -cwd
#$ -V
#$ -l h_vmem=4G
#$ -pe shm 12
#$ -l h_rt=5:59:00
#$ -l s_rt=5:59:00

## add modules & source specific conda environment
source activate CUT_n_RUN
module add ucsc_tools/2.7.2

## run script
python /srv/gsfs0/projects/snyder/chappell/scripts/CUTnRUN/PeakCalling/BAM2BDG_indiv.py \
--bamfile $BAMFILE \
--spikefile $SPIKEFILE \
--lengths_analysis $LENGTHS_ANALYSIS \
--lengths_image $LENGTHS_IMAGE \
--size_select_1 $SIZE_SELECT_1 \
--size_select_2 $SIZE_SELECT_2 \
--bedgraph $BEDGRAPH \
--chrom_sizes $CHROM_SIZES \
--big_wig $BIG_WIG \
--size_min_1 $SIZE_MIN_1 \
--size_max_1 $SIZE_MAX_1 \
--size_min_2 $SIZE_MIN_2 \
--size_max_2 $SIZE_MAX_2 \
--multiplying_factor $MULTIPLYING_FACTOR

## deactivate conda environment
source deactivate
EOF

## qsub then remove the tempscript
qsub $NAME.tempscript.sh
sleep 1
# rm $NAME.tempscript.sh

