#!/bin/bash

##### define your paths

# 1) a directory of fastq files
#fastq_files="/well/jknight/AbuDhabiRNA/P170245"

# 2) project directory
#mkdir -p ./mapping/2.mapping/mappingscripts
#mkdir -p ./mapping/2.mapping/fastq

# make links all fastq files in $fastq_files directory to $PROJECT_DIR/mapping/2.mapping/ ;
#cmd="ln -s $fastq_files/WTCHG_341122_20110*fastq.gz ./mapping/2.mapping/fastq/"
#echo $cmd


if [[ ! -e "mappingscripts" ]]; then

mkdir ./mappingscripts
fi


rst="./run.submit_jobs_for_mapping.sh"
echo "" > $rst

###################################################
#
# Add your sample data file here with a full path
 FILE="/well/jknight/AbuDhabiRNA/Katie/mapping/2.mapping/mapping.info.txt"
#
##################################################

while read -r line
do
	fastq_dir="$(echo "$line" | cut -f 1)"
	sample="$(echo "$line" | cut -f 1)"		# choose a right column for sample name
	echo $sample
	
	fastq="$(echo "$line" | cut -f 2)"
	fastq1="$(echo "$fastq" | sed -e 's/,/_1.fastq.gz,fastq\//g')"
	fastq2="$(echo "$fastq" | sed -e 's/,/_2.fastq.gz,fastq\//g')"
	fastq1=$fastq1"_1.fastq.gz"
	fastq2=$fastq2"_2.fastq.gz"
	fastq1="fastq/"$fastq1
	fastq2="fastq/"$fastq2

	echo $fastq1
	echo $fastq2


	exc_name="./mappingscripts/run.mapping.$sample.sh"
	cmd="sh 2.tophat2_mapping.sh $sample $fastq1 $fastq2"
	$cmd > $exc_name

	echo "qsub "$exc_name >> $rst
done < $FILE


### move all scripts to your working directory 
#mv ./mappingscripts/run.mapping* ./mapping/2.mapping/mappingscripts
#mv ./mappingscripts/run.submit_jobs_for_mapping.sh ./mapping/2.mapping




