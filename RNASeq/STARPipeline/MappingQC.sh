#\!/bin/bash
#$ -cwd -pe shmem 1 -V -N mappingQC
#$ -P jknight.prjc -q short.qc
#$ -o MappingQC.stdout
#$ -e MappingQC.sterr

module load R/3.2.2
cd $BASEDIR

Rscript --vanilla /well/jknight/RNASeqSTARPipeline/MappingQC.R
