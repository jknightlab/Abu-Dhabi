#!/bin/bash
# Make all scripts for the pipeline

##### move to your base directory containing the pipeline scripts and your 
# inital sample-fastq information key

cd $BASEDIR

# Make a sample key listing all the fastq files generated for each sample
module load R/3.2.2
Rscript --vanilla /well/jknight/RNASeqSTARPipeline/1.MakeSampleKey.R

# Make a script for each sample to align reads to the genome using STAR
sh /well/jknight/RNASeqSTARPipeline/2.MakeSTARMappingScripts.sh

# Make a script for each sample to filter out duplicates and output some mapping
# metrics using Picard
# sh /well/jknight/RNASeqSTARPipeline/MakePicardScripts.sh

# Make a script for each sample to count reads for each feature
# sh /well/jknight/RNASeqSTARPipeline/MakeCountingScripts.sh
