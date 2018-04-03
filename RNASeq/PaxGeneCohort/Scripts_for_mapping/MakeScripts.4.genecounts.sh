#!/bin/bash

source ~/.bashrc

# go to the job submission directory
working_dir="/well/jknight/AbuDhabiRNA/Katie/mapping/4.gene.counts/"
cd $working_dir


if [[ ! -e "countscripts" ]]; then

mkdir ./countscripts
fi


rst="run.submit_jobs_for_counting.sh"
echo "" > $rst

FILE="/well/jknight/AbuDhabiRNA/Katie/mapping/2.mapping/mapping.info.txt" 

while read -r line
do
        sample="$(echo "$line" | cut -f 1)"
        echo $sample


f_name="countscripts/run.job_"$sample".sh"

echo "#!/bin/bash" > $f_name
echo "#$ -cwd -pe shmem 2 -V -N HTSeqRead_"$sample >> $f_name
echo "#$ -P jknight.prjc -q short.qc" >> $f_name
echo "#$ -o /well/jknight/AbuDhabiRNA/Katie/mapping/4.gene.counts/"$sample".stdout" >> $f_name
echo "#$ -e /well/jknight/AbuDhabiRNA/Katie/mapping/4.gene.counts/"$sample".stderr" >> $f_name

echo "echo Start analysis at: \`date\`" >> $f_name

echo "export PATH=/apps/well/python/2.7/bin:/apps/well/python/2.7.6/bin:$PATH" >> $f_name
echo "export PYTHONPATH=/apps/well/python/2.7.6/bin:/well/jknight/software/rescomp/lib/python2.7/site-packages/:$PYTHONPATH" >> $f_name 
echo "export LD_LIBRARY_PATH=/apps/well/python/2.7.6/lib/:/apps/well/fftw/3.3.3-gcc4.7.2/lib:/apps/well/openblas/0.2.8-gcc4.7.2/lib:/apps/well/python/2.7.6/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" >> $f_name
echo "export PYTHONPATH=/apps/well/python/2.7.6/lib/python2.7/site-packages/:$PYTHONPATH" >> $f_name
echo "export PYTHONPATH=/well/jknight/Irina/Programs/HTSeq-0.6.1/build/lib.linux-x86_64-2.7:/well/jknight/software/rescomp/lib/python2.7/site-packages/${PYTHONPATH:+:$PYTHONPATH}" >> $f_name

cmd="/well/jknight/Irina/Programs/HTSeq-0.6.1/scripts/htseq-count --format=bam --order=pos --stranded=reverse /well/jknight/AbuDhabi/Katie/mapping/3.filtering/"$sample"/"$sample".nodup.NH1.properPairs.bam /well/jknight/reference/mapping/GRCh37/gencode.v19.chr_patch_hapl_scaff.annotation.gtf > /well/jknight/AbuDhabiRNA/Katie/mapping/4.gene_counts/"$sample".counts.txt"

echo $cmd >> $f_name 
echo "echo Finished analysis at: \`date\` " >> $f_name
echo "qsub "$f_name >> $rst

done < $FILE
