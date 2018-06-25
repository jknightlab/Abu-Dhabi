#!/bin/bash

echo "echo Start analysis at: \`date\`"

###################################################
#
# LOCATION OF REFERENCE AND TOOLS
#
###################################################

GFF=/well/jknight/AbuDhabiRNA/Katie/DiffExonUsage/Homo_sapiens.GRCh37.87.gff
SAMTOOLS="/apps/well/samtools/1.2/bin/samtools"

###################################################
#
# PROJECT DIRECTORY AND SAMPLE NAMES
#
###################################################

INPUT_DIR="/well/jknight/kate/RNAseq/P160233-REX/"
OUTPUT_DIR="/well/jknight/AbuDhabiRNA/Katie/DiffExonUsage/PilotData/"

SAMPLE_NAME=$1

###################################################
#
# PARAMETERS
#
###################################################

PAIRED="yes"
STRANDED="no"
FILETYPE="bam"

###################################################
#
# EXON COUNTING
#
###################################################


EXON_COUNTING(){
  
  echo "$SAMTOOLS sort -n ${INPUT_DIR}${SAMPLE_NAME}_hisat/accepted_hits.dedup.bam ${OUTPUT_DIR}${S1}.sorted"
  
  echo "$SAMTOOLS index ${OUTPUT_DIR}${S1}.sorted.bam"
  
  echo "python /well/jknight/AbuDhabiRNA/Katie/DiffExonUsage/dexseq_count.py $GFF -p $PAIRED -s $STRANDED -f $FILETYPE -r name ${OUTPUT_DIR}${S1}.sorted.bam ${OUTPUT_DIR}${S1}_exons.txt"

}


###################################################
#
# ENVIRONMENT SETUP
#
###################################################

ENVIRONMENT_SETUP(){
        echo "module load python/2.7.6"
        echo "export PYTHONPATH=/apps/well/python/2.7.6/bin:/well/jknight/software/rescomp/lib/python2.7/site-packages/:$PYTHONPATH"
        echo "echo"
}


###################################################
#
# WRITE A SCRIPT
#
###################################################

WRITE_A_SCRIPT(){

        echo "#\!/bin/bash"
        echo "#$ -cwd -pe shmem 1 -V -N exon."$1
        echo "#$ -P jknight.prjc -q short.qc"
        echo "#$ -o $OUTPUT_DIR$1.stdout.log"
        echo "#$ -e $OUTPUT_DIR$1.stderr.log"

        ENVIRONMENT_SETUP
	      EXON_COUNTING $1

        echo "echo Finished analysis at: \`date\`"

}

WRITE_A_SCRIPT $1
