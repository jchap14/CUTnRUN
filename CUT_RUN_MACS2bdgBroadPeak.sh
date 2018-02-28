#!/bin/bash
##### CUT&RUN: test for peaks using normalized human & fly bedgraphs

##### for x in `/bin/ls *K7ac.spikenorm.bdg` ; do bash CUT_RUN_MACS2bdgBroadPeak.sh $x; done

## load required modules
module add MACS2

## define variables
NAME=`basename $1 .spikenorm.bdg`

##### write a tempscript to be looped over
cat > $NAME.tempscript.sh << EOF
#!/bin/bash
#$ -N $NAME.broadPeak
#$ -j y
#$ -cwd
#$ -V
#$ -l h_vmem=2G
#$ -pe shm 12
#$ -l h_rt=5:59:00
#$ -l s_rt=5:59:00

##### run commands
## macs2 call peaks from bedgraph. cutoff is 1E-#
## help on https://github.com/taoliu/MACS/wiki/Advanced:-Call-peaks-using-MACS2-subcommands
## -l is minimum length of peak, better to set it as d value. Default is 200bp
## -g is max gap between significant points in a peak, better to set it as tag size. DEFAULT: 30

# -l 20 -g 30 gave 5-7x as many K27ac peaks relative to Zaugg EC set

macs2 bdgbroadcall -i $NAME.spikenorm.bdg --cutoff-peak 2 --cutoff-link 1 \
-l 20 -g 30 --lvl2-max-gap 100 --outdir . --o-prefix $NAME
EOF

## qsub then remove the tempscript
qsub $NAME.tempscript.sh 
sleep 1
rm $NAME.tempscript.sh