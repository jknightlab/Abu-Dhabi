#!/bin/bash

# USAGE:  sh.miseq_pipeline_testing.sh sample_name sample_A_1.fastq sample_A_2.fastq  > run_pipeline.sh
#        qsub run_pipeline.sh


echo "echo Start analysis at: \`date\`"

# Installed tools and references

SAMTOOLS="/apps/well/samtools/1.2/bin/samtools"
BEDTOOLS="/apps/well/bedtools/2.24.0-18-gb0bc5b7/bin/bedtools"
BEDGRAPH_TO_BIGWIG="/well/jknight/software/rescomp/bin/bedGraphToBigWig"
GENOME_COVERAGE_BED="/apps/well/bedtools/2.24.0-18-gb0bc5b7/bin/genomeCoverageBed"
BAMTOOLS="/apps/well/bamtools/2.3.0/bin/bamtools"
MARKDUP="/apps/well/picard-tools/1.111/MarkDuplicates.jar"



###################################################
#
# PROJECT DIRECTORY AND NAME OF DATA
#
###################################################

PROJECT_DIR="/well/jknight/AbuDhabiRNA/Katie/mapping/"		# you can define your project directory
INPUT_DIR=$PROJECT_DIR"2.mapping"
OUTPUT_DIR=$PROJECT_DIR"3.filtering"

if [[ ! -e $OUTPUT_DIR ]]; then
                echo "mkdir -p $OUTPUT_DIR"
fi

SAMPLE_NAME=$1



###################################################
#
# FILTERING Step1
#
###################################################


FILTERING_GENOME_MAPPING(){

	DIR_FILTER="$OUTPUT_DIR/$SAMPLE_NAME.pgf"
        if [[ ! -e $DIR_FILTER ]]; then
                echo "mkdir -p $DIR_FILTER"
        fi


	tmp="$SAMTOOLS sort $INPUT_DIR/$SAMPLE_NAME/accepted_hits.bam $DIR_FILTER/$SAMPLE_NAME.sorted"
	echo $tmp
	tmp="$SAMTOOLS index $DIR_FILTER/$SAMPLE_NAME.sorted.bam"
	echo $tmp

	tmp="java -Xmx8g -jar $MARKDUP I=$DIR_FILTER/$SAMPLE_NAME.sorted.bam O=$DIR_FILTER/$SAMPLE_NAME.nodup.bam REMOVE_DUPLICATES=true ASSUME_SORTED=true M=$DIR_FILTER/$SAMPLE_NAME.dup_metrix.txt"
	echo $tmp

	tmp="$BAMTOOLS filter -tag NH:1 -in $DIR_FILTER/$SAMPLE_NAME.nodup.bam -out $DIR_FILTER/$SAMPLE_NAME.nodup.NH1.bam"
        echo $tmp

	tmp="$SAMTOOLS view $DIR_FILTER/$SAMPLE_NAME.nodup.NH1.bam | cut -f 3,7 | sort | uniq -c > $DIR_FILTER/$SAMPLE_NAME.nodup.NH1.chrs.txt"
	echo $tmp

	tmp="$BAMTOOLS filter -isProperPair true -in $DIR_FILTER/$SAMPLE_NAME.nodup.NH1.bam -out $DIR_FILTER/$SAMPLE_NAME.nodup.NH1.properPairs.bam"
        echo $tmp
        tmp="$SAMTOOLS index $DIR_FILTER/$SAMPLE_NAME.nodup.NH1.properPairs.bam"
        echo $tmp

	tmp="$SAMTOOLS view $DIR_FILTER/$SAMPLE_NAME.nodup.NH1.properPairs.bam | cut -f 3,7 | sort | uniq -c > $DIR_FILTER/$SAMPLE_NAME.nodup.NH1.properPairs.chrs.txt"
        echo $tmp

}






###################################################
#
# ENVIRONMENT SETUP
#
###################################################

ENVIRONMENT_SETUP(){

        echo "module load R/3.1.3"
        echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/well/bamtools/2.3.0/lib/bamtools/"
        echo "sh /apps/well/profile.d/python/2.7.sh"
        echo "echo "
}



###################################################
#
# WRITE A SCRIPT
#
###################################################

WRITE_A_SCRIPT(){

        echo "#\!/bin/bash"
        echo "#$ -cwd -pe shmem 2 -V -N filter."$1
        echo "#$ -P jknight.prjc -q short.qc"
        echo "-o $OUTPUT_DIR$SAMPLE_NAME.pgf/log$1_stdout.log"
        echo "-e $OUTPUT_DIR$SAMPLE_NAME.pgf/log$1_stderr.log"

        ENVIRONMENT_SETUP
	FILTERING_GENOME_MAPPING $1

        echo "echo Finished analysis at: \`date\`"

}


WRITE_A_SCRIPT $1
