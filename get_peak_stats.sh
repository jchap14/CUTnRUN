#!/bin/bash
##### prep bamToBed output bed to format: Chr Start End Length

##### for x in `/bin/ls *.bed` ; do bash makeBed_chrStEndLength.sh $x; done

## define variables
name=`basename $1 .bed`

## write a tempscript to be looped over
cat > $name.tempscript.sh << EOF
#!/bin/bash
#$ -N $name.bedConvert
#$ -j y
#$ -cwd
#$ -V
#$ -l h_vmem=4G
#$ -pe shm 1
#$ -l h_rt=1:59:00
#$ -l s_rt=1:59:00

## run commands
mv $name.bed $name.oldBed
cat $name.oldBed | cut -f 1,2,3 | awk 'BEGIN { OFS = "\t" } { \$4 = \$3 - \$2 } 1' > $name.bed
EOF

## qsub then remove the tempscript
qsub $name.tempscript.sh 
sleep 1
rm $name.tempscript.sh
