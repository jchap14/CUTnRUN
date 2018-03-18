#!/bin/bash
# This script is how to submit the "BAM2BDG_indiv.py" script

##### submit for all sorted yeast bams in CWD (do in screen + QLOGIN for now))
## submit 
# cp /srv/gsfs0/projects/snyder/chappell/scripts/CUTnRUN/BAM2BDG_indiv.py .
# screen
# qlogin -l h_vmem=10G -l h_rt=24:00:00
# bash run_BAM2BDG.sh

## this doesn't work yet, for now use above to run on all BAMs at once
# for x in `/bin/ls *.nmSort.bam` ; do bash run_BAM2BDG.sh $x; done

## add modules
source activate CUT_n_RUN
module add ucsc_tools/2.7.2

## run script
python BAM2BDG.py

## deactivate conda environment
source deactivate

######## 

#!/bin/bash
# This is the submission script for "BAM2BDG_indiv.py"

# submit for a specific bedgraph
# bash run_BAM2BDG_indiv.sh $SAMPLE.nmSort.bam

##### submit for all sorted bedgraphs in CWD
# for x in `/bin/ls *.nmSort.bam` ; do bash run_BAM2BDG_indiv.sh $x; done

##### specify variables to pass to BAM2BDG_indiv.py
BAMFILE=$1
NAME=`basename $BAMFILE .nmSort.bam`
SPIKEFILE=`echo $NAME.yst.bam`
MERGE_CLOSE_PEAKS="True"
ends=False
lengths_analysis=True
lengths_image=True
size_select_1=True
size_select_2=True
bedgraph=True
chrom_sizes="hg19"
big_wig=True
size_min_1=20
size_max_1=120
size_min_2=150
size_max_2=710
multiplying_factor=10000

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
python /srv/gsfs0/projects/snyder/chappell/scripts/CUTnRUN/BAM2BDG_indiv.py \
                    --bedgraph $BEDGRAPH \
                    --threshold $THRESHOLD \
                    --min_length $MIN_LENGTH \
                    --inter_peak_distance $INTER_PEAK_DISTANCE \
                    --merge_close_peaks $MERGE_CLOSE_PEAKS \
                    --max_length $MAX_LENGTH \
                    --generate_ID $GENERATE_ID

python /srv/gsfs0/projects/snyder/chappell/scripts/CUTnRUN/PeakCalling/BAM2BDG.py \

## deactivate conda environment
source deactivate
EOF

## qsub then remove the tempscript
qsub $NAME.tempscript.sh
sleep 1
rm $NAME.tempscript.sh

