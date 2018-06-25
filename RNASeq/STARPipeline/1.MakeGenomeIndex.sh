#!/bin/bash
#$ -cwd -pe shmem 4 -V -N index
#$ -P jknight.prjc -q short.qc

echo Start analysis at: `date`

export PATH=$PATH:/apps/well/star/20151002/bin/Linux_x86_64

# https://leonjessen.wordpress.com/2014/12/01/how-do-i-create-star-indices-using-the-newest-grch38-version-of-the-human-genome/

# mkdir -p data/GRCh38/sequence
# cd data/GRCh38/sequence/
# wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz
# wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{MT,X,Y}.fa.gz
# gunzip -c Homo_sapiens.GRCh38.dna.chromosome.* > GRCh38_r92.all.fa
# cd ../../../
# 
# mkdir -p data/GRCh38/annotation
# cd data/GRCh38/annotation/
# wget ftp://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.gtf.gz
# gunzip Homo_sapiens.GRCh38.92.gtf.gz
# cd ../../../

mkdir -p /well/jknight/RNASeqSTARPipeline/data/GRCh38/star_indices_overhang74
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /well/jknight/RNASeqSTARPipeline/data/GRCh38/star_indices_overhang74/ --genomeFastaFiles data/GRCh38/sequence/GRCh38_r92.all.fa --sjdbGTFfile data/GRCh38/annotation/Homo_sapiens.GRCh38.92.gtf --sjdbOverhang 74

echo Finished analysis at: `date`

/apps/well/star/20151002/bin/Linux_x86_64/STAR --runThreadN 1 --runMode genomeGenerate --genomeDir /well/jknight/RNASeqSTARPipeline/data/GRCh38/star_indices_overhang74 --genomeFastaFiles data/GRCh38/sequence/GRCh38_r92.all.fa --sjdbGTFfile data/GRCh38/annotation/Homo_sapiens.GRCh38.92.gtf --sjdbOverhang 74