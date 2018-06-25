#!/bin/bash

# http://htseq.readthedocs.io/en/release_0.9.1/count.html?highlight=htseq-count

FILE=$BASEDIR"/mapping.info.txt"

# Make a directory for the counting step and move into it
if [[ ! -e "Counts" ]]; then

mkdir ./Counts
fi

cd Counts

OUTPUT_DIR=$BASEDIR"/Counts"

# Make a directory for the filtering scripts
if [[ ! -e "Countscripts" ]]; then

mkdir ./Countscripts
fi

# Start a script that will submit all the per-sample counting scripts
rst="./submit_jobs_for_counting.sh"
echo "" > $rst

# read in the sample key
# for each sample input the bam file to the counting script
while read -r line
do
	sample="$(echo "$line" | cut -f 1)"
	echo $sample
	
	exc_name="$OUTPUT_DIR/Countscripts/count."$sample".sh"
	
	echo "#!/bin/bash" >> $exc_name
  echo "#$ -cwd -pe shmem 1 -V -N count."$sample >> $exc_name
  echo "#$ -P jknight.prjc -q short.qc" >> $exc_name
  echo "#$ -o $OUTPUT_DIR/$sample.stdout" >> $exc_name
  echo "#$ -e $OUTPUT_DIR/$sample.stderr" >> $exc_name

  echo "echo Start analysis at: \`date\`" >> $exc_name

  echo "export PATH=/apps/well/python/2.7/bin:/apps/well/python/2.7.6/bin:$PATH" >> $exc_name
  echo "export PYTHONPATH=/apps/well/python/2.7.6/bin:/well/jknight/software/rescomp/lib/python2.7/site-packages/:$PYTHONPATH" >> $exc_name
  echo "module load python/2.7.6" >> $exc_name
  echo "export LD_LIBRARY_PATH=/apps/well/python/2.7.6/lib/:/apps/well/fftw/3.3.3-gcc4.7.2/lib:/apps/well/openblas/0.2.8-gcc4.7.2/lib:/apps/well/python/2.7.6/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" >> $exc_name
  echo "export PYTHONPATH=/apps/well/python/2.7.6/lib/python2.7/site-packages/:$PYTHONPATH" >> $exc_name
  echo "export PYTHONPATH=/well/jknight/Irina/Programs/HTSeq-0.6.1/build/lib.linux-x86_64-2.7:/well/jknight/software/rescomp/lib/python2.7/site-packages/${PYTHONPATH:+:$PYTHONPATH}" >> $exc_name

  cmd="/well/jknight/RNASeqSTARPipeline/htseq-count --type=exon --order=pos --mode=union --stranded=no --format=bam --idattr=gene_name $BASEDIR/Filtering/"$sample"/"$sample".nodup.NH1.properPairs.bam $GTF > $OUTPUT_DIR/$sample.counts.txt"
  echo $cmd >> $exc_name 
  echo "echo Finished analysis at: \`date\` " >> $exc_name
  
  echo "qsub -hold_jid filter."$sample $exc_name >> $rst

done < $FILE

cd ..






