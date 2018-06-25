#!/bin/bash

# Run pipeline

module load R/3.2.2

CONFIG=$1

echo "Using configuration file" $CONFIG

echo "Setting up pipeline using your configuration file" $CONFIG

source $CONFIG

export GTF=/well/jknight/RNASeqSTARPipeline/data/GRCh38/annotation/Homo_sapiens.GRCh38.92.gtf

export GENOME=/well/jknight/RNASeqSTARPipeline/data/GRCh38/star_indices_overhang74

if [ "$READLENGTH" -ne "75" ]
then
  echo "The pipeline will be run using the generic value of 100 for the genome index. According to the STAR manual, this will work well but is not optimised for your read length."
  export GENOME=/well/jknight/RNASeqSTARPipeline/data/GRCh38/star_indices_overhang100
fi

##### move to your base directory containing the pipeline scripts and your 
# inital sample-fastq information key

cd $BASEDIR

echo "Making sample key"

# Make a sample key listing all the fastq files generated for each sample
module load R/3.2.2
Rscript --vanilla /well/jknight/RNASeqSTARPipeline/1.MakeSampleKey.R $KEY

echo "Making scripts for each sample"

# Make a script for each sample to align reads to the genome using STAR
sh /well/jknight/RNASeqSTARPipeline/2.MakeSTARMappingScripts.sh

# Make a script for each sample to filter out duplicates and output some mapping
# metrics using Picard
sh /well/jknight/RNASeqSTARPipeline/3.MakeFilteringScripts.sh

# Make a script for each sample to count reads for each feature
sh /well/jknight/RNASeqSTARPipeline/4.MakeCountScripts.sh

echo "Submitting jobs for alignment"

sh $BASEDIR/Mapping/submit_jobs_for_mapping.sh
rm $BASEDIR/Mapping/submit_jobs_for_mapping.sh

echo "Mapping QC"  # hold on mapping

qsub -hold_jid "map.*" /well/jknight/RNASeqSTARPipeline/MappingQC.sh

echo "Submitting jobs for filtering"

sh $BASEDIR/Filtering/submit_jobs_for_filtering.sh
rm $BASEDIR/Filtering/submit_jobs_for_filtering.sh

echo "Submitting jobs for counting"

sh $BASEDIR/Counts/submit_jobs_for_counting.sh
rm $BASEDIR/Counts/submit_jobs_for_counting.sh

echo "Merge counts"

qsub -hold_jid "count.*" /well/jknight/RNASeqSTARPipeline/5.MergeCounts.sh
