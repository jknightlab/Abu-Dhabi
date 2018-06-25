#\!/bin/bash
#$ -cwd -pe shmem 1 -V -N gff
#$ -P jknight.prjc -q short.qc
#$ -o gff_stdout.log
#$ -e gff_stderr.log

# Substitute the appropriate annotation file
# Here is the one used in Wan's RNA-Seq mapping

python /well/jknight/AbuDhabiRNA/Katie/DiffExonUsage/dexseq_prepare_annotation.py \
/well/jknight/AbuDhabiRNA/Katie/mapping/2.mapping/symlink/gencode.v19.chr_patch_hapl_scaff.annotation.gtf \
gencode.v19.chr_patch_hapl_scaff.annotation.gff