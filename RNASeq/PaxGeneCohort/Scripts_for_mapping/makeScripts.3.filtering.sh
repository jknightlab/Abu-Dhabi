#!/bin/bash


if [[ ! -e "filteringscripts" ]]; then

mkdir ./filteringscripts
fi


rst="./run.submit_jobs_for_filtering.sh"
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
	sample="$(echo "$line" | cut -f 1)"		# choose the right column for sample name
	echo $sample
	

	exc_name="./filteringscripts/run.filtering.$sample.sh"
	cmd="sh 3.2.filtering.genome.sh $sample"
	$cmd > $exc_name

	echo "qsub "$exc_name >> $rst
done < $FILE

