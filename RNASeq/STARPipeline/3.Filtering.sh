#!/bin/bash

###################################################
#
#  SETUP PATHS TO TOOLS AND REFERENCE
#
###################################################

export PATH=$PATH:/apps/well/samtools/1.2/bin
SAMTOOLS="/apps/well/samtools/1.2/bin/samtools"
export PATH=$PATH:/apps/well/bedtools/2.24.0-18-gb0bc5b7/bin/
BEDTOOLS="/apps/well/bedtools/2.24.0-18-gb0bc5b7/bin/bedtools"
BEDGRAPH_TO_BIGWIG="/well/jknight/software/rescomp/bin/bedGraphToBigWig"
GENOME_COVERAGE_BED="/apps/well/bedtools/2.24.0-18-gb0bc5b7/bin/genomeCoverageBed"
BAMTOOLS="/apps/well/bamtools/2.3.0/bin/bamtools"
export PATH=$PATH:/apps/well/bamtools/2.3.0/bin
export PATH=$PATH:/apps/well/bamtools/2.3.0/lib/bamtools
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/well/bamtools/2.3.0/lib/bamtools/
MARKDUP="/apps/well/picard-tools/1.111/MarkDuplicates.jar"


###################################################
#
# PROJECT DIRECTORY AND INPUT DATA
#
###################################################

INPUT_DIR=$BASEDIR"/Mapping/"

OUTPUT_DIR=$BASEDIR"/Filtering/"

SAMPLE_NAME=$1


###################################################
#
# FILTERING
#
###################################################

DIR_SAMPLE_NAME=$OUTPUT_DIR$SAMPLE_NAME

FILTERING(){
  
  echo
	echo "echo Start filtering at: \`date\`"
  echo

	tmp="$SAMTOOLS sort $INPUT_DIR$SAMPLE_NAME/${SAMPLE_NAME}Aligned.out.bam $DIR_SAMPLE_NAME/$SAMPLE_NAME.sorted"
	echo $tmp
	
	# tmp="rm $INPUT_DIR/$SAMPLE_NAME/$SAMPLE_NAMEAligned.out.bam"
	# echo $tmp

	tmp="$SAMTOOLS index $DIR_SAMPLE_NAME/$SAMPLE_NAME.sorted.bam"
	echo $tmp
	
	tmp="java -Xmx8g -jar $MARKDUP I=$DIR_SAMPLE_NAME/$SAMPLE_NAME.sorted.bam O=$DIR_SAMPLE_NAME/$SAMPLE_NAME.nodup.bam REMOVE_DUPLICATES=true ASSUME_SORTED=true M=$DIR_SAMPLE_NAME/$SAMPLE_NAME.dup_metrix.txt"
	echo $tmp
	
	tmp="rm $DIR_SAMPLE_NAME/$SAMPLE_NAME.sorted.bam"
	echo $tmp

	tmp="$BAMTOOLS filter -tag NH:1 -in $DIR_SAMPLE_NAME/$SAMPLE_NAME.nodup.bam -out $DIR_SAMPLE_NAME/$SAMPLE_NAME.nodup.NH1.bam"
  echo $tmp
  
	tmp="rm $DIR_SAMPLE_NAME/$SAMPLE_NAME.nodup.bam"
	echo $tmp

	tmp="$SAMTOOLS view $DIR_SAMPLE_NAME/$SAMPLE_NAME.nodup.NH1.bam | cut -f 3,7 | sort | uniq -c > $DIR_SAMPLE_NAME/$SAMPLE_NAME.nodup.NH1.chrs.txt"
	echo $tmp

	tmp="$BAMTOOLS filter -isProperPair true -in $DIR_SAMPLE_NAME/$SAMPLE_NAME.nodup.NH1.bam -out $DIR_SAMPLE_NAME/$SAMPLE_NAME.nodup.NH1.properPairs.bam"
  echo $tmp

	tmp="rm $DIR_SAMPLE_NAME/$SAMPLE_NAME.nodup.NH1.bam"
	echo $tmp
  
  tmp="$SAMTOOLS index $DIR_SAMPLE_NAME/$SAMPLE_NAME.nodup.NH1.properPairs.bam"
  echo $tmp

	tmp="$SAMTOOLS view $DIR_SAMPLE_NAME/$SAMPLE_NAME.nodup.NH1.properPairs.bam | cut -f 3,7 | sort | uniq -c > $DIR_SAMPLE_NAME/$SAMPLE_NAME.nodup.NH1.properPairs.chrs.txt"
  echo $tmp
        
	echo "echo FILTERING DONE."
}

###################################################
#
# ENVIRONMENT SETUP
#
###################################################

ENVIRONMENT_SETUP(){

        echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/well/bamtools/2.3.0/lib/bamtools/"
        echo "echo "

}


###################################################
#
# WRITE A SCRIPT
#
###################################################

WRITE_A_SCRIPT(){

        echo "#\!/bin/bash"
        echo "#$ -cwd -pe shmem 1 -V -N filter."$1
        echo "#$ -P jknight.prjc -q short.qc"
        echo "#$ -o $DIR_SAMPLE_NAME/$1.stdout"
        echo "#$ -e $DIR_SAMPLE_NAME/$1.sterr"
	
        
	ENVIRONMENT_SETUP
	FILTERING $1
	
        echo "echo Finished analysis at: \`date\`"

}

WRITE_A_SCRIPT $1
