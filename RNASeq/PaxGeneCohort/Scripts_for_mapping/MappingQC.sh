#\!/bin/bash
#$ -cwd -pe shmem 2 -V -N rnaseqc
#$ -P jknight.prjc -q short.qc
#$ -o qc_stdout.log
#$ -e qc_stderr.log

echo Start analysis at: `date`

source ~/.bashrc
echo 
echo =========================================

java -jar /apps/well/rna-seqc/1.1.8/RNA-SeQC_v1.1.8.jar \
-o /well/jknight/AbuDhabiRNA/Katie/QC/RNASeQC/ \
-r /well/jknight/AbuDhabiRNA/Katie/mapping/2.mapping/symlink/Gencode19.GRCh37.genome.PGF.decoy.ERCC.fa \
-s /well/jknight/AbuDhabiRNA/Katie/QC/bamfiles.txt \
-t /well/jknight/AbuDhabiRNA/Katie/mapping/2.mapping/symlink/gencode.v19.chr_patch_hapl_scaff.annotation.gtf

echo Finished analysis at: `date`