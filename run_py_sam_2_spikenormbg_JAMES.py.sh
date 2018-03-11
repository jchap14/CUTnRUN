#!/bin/bash
# This script is how to submit the "run_py_sam_2_spikenormbg_JAMES.py.sh" script

##### submit for all sorted yeast bams in CWD (do in screen QLOGIN for now))
## submit 
# cp /srv/gsfs0/projects/snyder/chappell/scripts/CUT_and_RUN_scripts/py_sam_2_spikenormbg_JAMES.py .
# bash run_py_sam_2_spikenormbg_JAMES.py.sh

## this doesn't work yet, for now use above to run on all BAMs at once
# for x in `/bin/ls *.yst.sorted.bam` ; do bash run_py_sam_2_spikenormbg_JAMES.py.sh $x; done

## add modules
source activate CUT_n_RUN
module add ucsc_tools/2.7.2

## run script
python py_sam_2_spikenormbg_JAMES.py

## deactivate conda environment
source deactivate

## change all .bg files to .bdg

