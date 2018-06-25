#!/bin/bash

FILE=$BASEDIR"/mapping.info.txt"

# Make a directory for the filtering step and move into it
if [[ ! -e "Filtering" ]]; then

mkdir ./Filtering
fi

cd Filtering

OUTPUT_DIR=$BASEDIR"/Filtering/"

# Make a directory for the filtering scripts
if [[ ! -e "Filteringscripts" ]]; then

mkdir ./Filteringscripts
fi

# Start a script that will submit all the per-sample filtering scripts
rst="./submit_jobs_for_filtering.sh"
echo "" > $rst

# read in the sample key
# for each sample input the bam file to the filtering script
while read -r line
do
	sample="$(echo "$line" | cut -f 1)"
	echo $sample
	
	DIR_SAMPLE_NAME=$OUTPUT_DIR$sample
	if [[ ! -e $DIR_SAMPLE_NAME ]]; then
                mkdir -p $DIR_SAMPLE_NAME
  fi

	exc_name=$OUTPUT_DIR"Filteringscripts/filtering."$sample.sh
	cmd="sh /well/jknight/RNASeqSTARPipeline/3.Filtering.sh $sample"
	$cmd > $exc_name

	echo "qsub -hold_jid map."$sample $exc_name >> $rst
done < $FILE

cd ..
