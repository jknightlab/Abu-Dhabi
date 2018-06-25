#!/bin/bash

###################################################
#
#  SETUP PATHS TO TOOLS AND REFERENCE
#
###################################################

module load python/2.7.6
export PYTHONPATH=/apps/well/python/2.7.6/bin:/well/jknight/software/rescomp/lib/python2.7/site-packages/:$PYTHONPATH
export PATH=$PATH:/apps/well/bowtie2/2.2.5
export PATH=$PATH:/apps/well/samtools/1.2/bin
export PATH=$PATH:/apps/well/bedtools/2.24.0-18-gb0bc5b7/bin/
export PATH=$PATH:/apps/well/bamtools/2.3.0/bin
export PATH=$PATH:/apps/well/bamtools/2.3.0/lib/bamtools
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/well/bamtools/2.3.0/lib/bamtools/
export PATH=$PATH:/apps/well/star/20151002/bin/Linux_x86_64
export STAR=/well/jknight/RNASeqSTARPipeline/STAR-2.6.0a/bin/Linux_x86_64/STAR

###################################################
#
# PROJECT DIRECTORY AND INPUT DATA
#
###################################################

OUTPUT_DIR=$BASEDIR"/Mapping/"

SAMPLE_NAME=$1

INPUT_READ1=$2
INPUT_READ2=$3

###################################################
#
# STAR MAPPING
#
###################################################

DIR_SAMPLE_NAME=$OUTPUT_DIR$SAMPLE_NAME

STAR_MAPPING(){
        echo
	      echo "echo Start mapping at: \`date\`"
        echo

	tmp="$STAR --genomeDir $GENOME \
	          --runThreadN 4 \
	          --readFilesIn $INPUT_READ1 $INPUT_READ2 \
	        --readFilesCommand gunzip -c \
		--outFileNamePrefix $DIR_SAMPLE_NAME"/"$SAMPLE_NAME \
	          --outSAMtype BAM Unsorted \
	          --outFilterMismatchNoverLmax 0.067 \
	          --outSAMprimaryFlag OneBestScore \
	          --outReadsUnmapped Fastx 2> $DIR_SAMPLE_NAME/stderr.$SAMPLE_NAME.STAR.txt"
	echo $tmp

	echo "echo MAPPING DONE."	
}

###################################################
#
# ENVIRONMENT SETUP
#
###################################################

ENVIRONMENT_SETUP(){

        echo "module load R/3.2.2"
        echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/well/bamtools/2.3.0/lib/bamtools/"
        echo "sh /apps/well/profile.d/python/2.7.sh"
	      echo "export PATH=$PATH:/apps/well/star/20151002/bin/Linux_x86_64"
        echo "echo "
}


###################################################
#
# WRITE A SCRIPT
#
###################################################

WRITE_A_SCRIPT(){

        echo "#\!/bin/bash"
        echo "#$ -cwd -pe shmem 4 -V -N map."$1
        echo "#$ -P jknight.prjc -q short.qc"
        echo "#$ -o $DIR_SAMPLE_NAME/$1.stdout"
        echo "#$ -e $DIR_SAMPLE_NAME/$1.sterr"
	
        
#	ENVIRONMENT_SETUP
	STAR_MAPPING $1 $2 $3
	
        echo "echo Finished analysis at: \`date\`"

}

WRITE_A_SCRIPT $1 $2 $3 $4
