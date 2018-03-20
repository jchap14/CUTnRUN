#!/bin/bash
# This is the submission script for "py_peak_calling_JAMES.py"

# submit for a specific bedgraph
# bash run_py_peak_calling_JAMES.py.sh $SAMPLE.bdg

##### submit for all sorted bedgraphs in CWD
# for x in `/bin/ls *.nmSort.bdg` ; do bash run_py_peak_calling_JAMES.py.sh $x; done

##### specify variables to pass to py_peak_calling_James.py
BEDGRAPH=$1
## THRESHOLD is not an FDR, but score based?
THRESHOLD=10
MIN_LENGTH=50
INTER_PEAK_DISTANCE=100
MERGE_CLOSE_PEAKS="True"
MAX_LENGTH=1000
GENERATE_ID="True"

NAME=`basename $BEDGRAPH .nmSort.bdg`

## create tempscript
cat > $NAME.tempscript.sh << EOF
#!/bin/bash
#$ -N $NAME.$THRESHOLD.PeakCall
#$ -S /bin/bash
#$ -j y
#$ -cwd
#$ -V
#$ -l h_vmem=4G
#$ -pe shm 1
#$ -l h_rt=0:59:00
#$ -l s_rt=0:59:00

## add modules & source specific conda environment
source activate CUT_n_RUN

## run script
python /srv/gsfs0/projects/snyder/chappell/JR/SLO/CUT_n_RUN/original_set/PeakCalling/py_peak_calling_JAMES.py \
                    --bedgraph $BEDGRAPH \
                    --threshold $THRESHOLD \
                    --min_length $MIN_LENGTH \
                    --inter_peak_distance $INTER_PEAK_DISTANCE \
                    --merge_close_peaks $MERGE_CLOSE_PEAKS \
                    --max_length $MAX_LENGTH \
                    --generate_ID $GENERATE_ID

## deactivate conda environment
source deactivate
EOF

## qsub then remove the tempscript
qsub $NAME.tempscript.sh
sleep 1
rm $NAME.tempscript.sh

