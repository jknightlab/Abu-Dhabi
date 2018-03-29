#!/bin/bash


# USAGE: 2.tophat2_mapping.sh  sample_name project_dir sample_A_1.fastq sample_A_2.fastq  > run_pipeline.sh
# USAGE: pipeline_tophat2_mapping.sh  sample_name sample_A_1.fastq,sample_B_1.fastq sample_A_2.fastq,sample_B_2.fastq  > run_pipeline.sh
#        qsub run_pipeline.sh


echo "echo Start analysis at: \`date\`"


###################################################
#
#  SETUP PATHS
#
###################################################

module load R/3.1.3
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/well/bamtools/2.3.0/lib/bamtools/
sh /apps/well/profile.d/python/2.7.sh

# Installed tools and references

# BOWTIE="/apps/well/bowtie2/2.2.5" # make sure this version of bowtie for compatibility with tophat
SEQTK="/well/jknight/Wan/tools/cluster3/bin/seqtk"
SAMTOOLS="/apps/well/samtools/1.2/bin/samtools"
BEDTOOLS="/apps/well/bedtools/2.24.0-18-gb0bc5b7/bin/bedtools"
BAMTOOLS="/apps/well/bamtools/2.3.0/bin/bamtools"

TOPHAT2=/apps/well/tophat/2.0.14/bin/tophat2

# GRCh37
GTF=/well/jknight/AbuDhabiRNA/Katie/mapping/2.mapping/symlink/gencode.v19.chr_patch_hapl_scaff.annotation
BOWTIE2_GENOME_INDEX=/well/jknight/AbuDhabiRNA/Katie/mapping/2.mapping/symlink/Gencode19.GRCh37.genome.PGF.decoy.ERCC

###################################################
#
# TOPHAT2 PARAMETERS
#
###################################################

nthreads=4			# default: 1 
mate_inner_dist=190		# default: 50
mate_std_dev=50		# default: 20
read_mismatches=4		# default: 2
read_edit_dist=4		# default: 2
b2_N=1				# default: 0 multiseed alignment number of mismatches
library_type="fr-unstranded"	# default is unstranded (fr-unstranded)

read_gap_length=2		# default: 2
segment_mismatches=2		# default: 2

###################################################
#
# PROJECT DIRECTORY AND NAME OF DATA
#
###################################################

PROJECT_DIR="/well/jknight/AbuDhabiRNA/Katie/"	# you can define your project directory
SAMPLE_NAME=$1

OUTPUT_DIR=$PROJECT_DIR"mapping/2.mapping/"

if [[ ! -e $OUTPUT_DIR ]]; then
                mkdir -p $OUTPUT_DIR
fi

SAMPLE_NAME=$1

###################################################
#
# INPUT DATA: fastq files
#
###################################################

INPUT_READ1=$2
INPUT_READ2=$3

###################################################
#
# TOPHAT2 MAPPING
#
###################################################

TOPHAT2_MAPPING(){

	echo "echo ========================================="
	echo "echo TOPHAT2 MAPPING: \`date\`"

	DIR_SAMPLE_NAME=$OUTPUT_DIR$SAMPLE_NAME
	if [[ ! -e $DIR_SAMPLE_NAME ]]; then
		echo "mkdir -p $DIR_SAMPLE_NAME"
	fi

	#tmp="cd $OUTPUT_DIR"
	#echo $tmp


	tmp="$TOPHAT2 -p $nthreads \
			--b2-very-sensitive \
			--mate-inner-dist $mate_inner_dist \
			--mate-std-dev $mate_std_dev \
			--read-mismatches $read_mismatches \
			--read-edit-dist $read_edit_dist \
			--b2-N 1 \
			--library-type $library_type \
			--read-gap-length $read_gap_length \
			--segment-mismatches $segment_mismatches \
			-G $GTF.gtf \
			--transcriptome-index $GTF \
			-o $SAMPLE_NAME $BOWTIE2_GENOME_INDEX $2 $3 \
			2> stderr.$SAMPLE_NAME.tophat2.txt"
	echo $tmp


	echo "Manipulations on sam file:  indexing: \`date\`"

	tmp="$SAMTOOLS index $OUTPUT_DIR$SAMPLE_NAME/accepted_hits.bam \
		2> stderr.$SAMPLE_NAME.samtools_index.txt"
	echo $tmp

	echo "echo MAPPING DONE."	
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
	      echo "source ~/.bashrc"
        echo "echo "
}




###################################################
#
# WRITE A SCRIPT
#
###################################################

WRITE_A_SCRIPT(){

        echo "#\!/bin/bash"
        echo "#$ -cwd -pe shmem 4 -V -N log."$1
        echo "#$ -P jknight.prjc -q long.qc"
        echo "#$ -o $1/$1.stdout"
        echo "#$ -e $1/$1.sterr"

        ENVIRONMENT_SETUP
	      TOPHAT2_MAPPING $1 $2 $3
	
        echo "echo Finished analysis at: \`date\`"

}


WRITE_A_SCRIPT $1 $2 $3 $4
