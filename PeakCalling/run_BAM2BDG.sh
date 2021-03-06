#!/bin/bash
# This script is how to submit the "run_py_sam_2_spikenormbg_JAMES.py.sh" script

##### submit for all sorted yeast bams in CWD (do in screen + QLOGIN for now))
## submit 
# cp /srv/gsfs0/projects/snyder/chappell/scripts/CUTnRUN/BAM2BDG.py .
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

## change all .bg files to .bdg

