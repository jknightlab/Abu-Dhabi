#\!/bin/bash
#$ -cwd -pe shmem 1 -V -N mergecounts
#$ -P jknight.prjc -q short.qc
#$ -o mergecounts.stdout
#$ -e mergecounts.sterr

module load R/3.2.2
cd $BASEDIR

Rscript --vanilla /well/jknight/RNASeqSTARPipeline/5.MergeCounts.R

rm -R ./Counts
