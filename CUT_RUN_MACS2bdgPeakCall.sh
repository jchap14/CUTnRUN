#!/bin/bash
##### CUT&RUN: test for peaks using yeast normalized human bedgraphs

##### for x in `/bin/ls *.nmSort.bdg` ; do bash CUT_RUN_MACS2bdgPeakCall.sh {CUTOFF} $x; done

## load required modules
module add MACS2

## define variables
CUTOFF=$1 #cutoff is 1E^-number
BEDGRAPH=$2
NAME=`basename $1 .nmSort.bdg`

##### set macs2 bdgpeakcall options -l & -g
## -l is minimum length of peak, better to set it as d value. Default is 200bp
## -g is max gap between significant points in a peak, better to set it as tag size. DEFAULT: 30

# -l 20 -g 30 gave 5-7x as many K27ac peaks relative to Zaugg EC set
L_OPTION=50
G_OPTION=100
OUTPREFIX=`echo $NAME."FDRe"$CUTOFF`

##### write a tempscript to be looped over
cat > $NAME.tempscript.sh << EOF
#!/bin/bash
#$ -N $OUTPREFIX
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

macs2 bdgpeakcall -i $BEDGRAPH --cutoff $CUTOFF -l $L_OPTION -g $G_OPTION \
--no-trackline --outdir . --o-prefix $OUTPREFIX
EOF

## qsub then remove the tempscript
qsub $NAME.tempscript.sh 
sleep 1
rm $NAME.tempscript.sh