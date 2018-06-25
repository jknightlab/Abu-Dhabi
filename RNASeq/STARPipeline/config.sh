#!/bin/bash

# This is your configuration file
# It is the only file you need to edit
# Put it in your project directory
# When you run the pipeline use the path to this file as an argument
# i.e. sh /well/jknight/RNASeqSTARPipeline/RunPipeline.sh /well/jknight/AbuDhabiRNA/Katie/STAR

###################################
#
# Paths to your files
#
###################################

# Your project directory
export BASEDIR=/well/jknight/AbuDhabiRNA/Katie/STAR
echo "Base directory: "$BASEDIR

# Pathway to your sample information file
export FILE=/well/jknight/AbuDhabiRNA/Katie/STAR/sample_fastq_key.txt
echo "Sample information file: "$FILE

# Pathway to directory containing your fastq files (or symlinks to them)
export FASTQ=/well/jknight/AbuDhabiRNA/Katie/mapping/2.mapping/fastq
echo "Fastq directory: "$FASTQ

# Path to STAR
export STAR=/well/jknight/RNASeqSTARPipeline/STAR-2.6.0a/bin/Linux_x86_64/STAR
echo "STAR: "$STAR

###################################
#
# STAR mapping parameters
#
###################################

export READLENGTH=75
echo "Read length: "$READLENGTH