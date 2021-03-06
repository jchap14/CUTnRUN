#!/bin/bash

##### Submission script for BDS controlled ChIPseq pipeline
## for x in `/bin/ls *.trim.R1.fq.gz` ; do bash Kundaje_ChIPseq_CnR.sh $x; done

##### INPUTs required
## type of factor for ChIP (-type in python script)
CHIPTYPE="histone"
# CHIPTYPE="TF"

##### input files to pass to chipseq.py (FQs should be trimmed)
## calc replicate 1 FQs (format == $NAME.repB)
FQ_R1=$1
NAME=`basename $FQ_R1 .trim.R1.fq.gz`
FQ_R2=`echo $NAME.trim.R2.fq.gz`

## create tempscript
cat > $NAME.tempscript.sh << EOF
#!/bin/bash

## conda environment
## don't source conda envs, they are auto-activated by the bds pipeline

##### run script
## -type can be histone or TF
## forcing macs2 for peak calling
## maybe xcor can be turned off, but then need to set -extsize_macs2 [EXTSIZE]

python /srv/gsfs0/projects/snyder/chappell/TF_chipseq_pipeline/chipseq.py \
-type $CHIPTYPE --screen $NAME -system slurm -q_for_slurm_account -q mpsnyder \
-pe -species hg19 -nth 12 -peak_caller macs2 \
-fastq1_1 $FQ_R1 -fastq1_2 $FQ_R2 \
-out_dir $NAME -mem_dedup 80G -mem_bwa 50G -mem_spp 20G -mem_macs2 25G

## deactivate conda environment
EOF

## qsub then remove the tempscript
bash $NAME.tempscript.sh
sleep 1
rm $NAME.tempscript.sh
