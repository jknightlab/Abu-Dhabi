#!/bin/bash


if [[ ! -e "exoncountscripts" ]]; then

mkdir ./exoncountscripts
fi


rst="./run.submit_jobs_for_exon_counting.sh"
echo "" > $rst



###################################################
#
# Add your sample data file here with a full path
 FILE="/well/jknight/kate/RNAseq/P160233-REX/README.txt"
#
##################################################


while read -r line
do
	fastq_dir="$(echo "$line" | cut -f 1)"
	sample="$(echo "$line" | cut -f 1)"		# choose the right column for sample name
	echo $sample
	

	exc_name="./exoncountscripts/run.exon.counts.$sample.sh"
	cmd="sh exon_counting.sh $sample"
	$cmd > $exc_name

	echo "qsub "$exc_name >> $rst
done < $FILE

