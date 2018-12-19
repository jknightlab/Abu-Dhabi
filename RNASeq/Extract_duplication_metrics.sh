#!/bin/bash

### Pull out Picard duplication metrics

FILE="/well/jknight/Core_RNASeq_Input_Trial/mapping.info.txt"

while read -r line
do

	sample="$(echo "$line" | cut -f 1)"
	echo $sample
	DIR_SAMPLE_NAME=/well/jknight/Core_RNASeq_Input_Trial/Filtering/$sample
  cd $DIR_SAMPLE_NAME
  grep Unknown $sample.dup_metrix.txt > $sample.dup.counts.txt

done < $FILE

cd ..





R
library(data.table)
log.paths <- list.files(path=".", pattern="*dup.counts.txt", full.names=TRUE, recursive=TRUE)
sample.names <- list.files(path=".", pattern="*dup.counts.txt", full.names=FALSE, recursive=TRUE)
sample.names <- unlist(strsplit(sample.names, "/"))[seq(1, 2*length(sample.names), 2)]
logs <- lapply(log.paths, read.delim, header=FALSE)
logs <- lapply(logs, as.data.frame)
dup.logs <- rbindlist(logs)
rownames(dup.logs) <- sample.names
colnames(dup.logs) <- c("Library", "Unpaired_Reads_Examined", "Read_Pairs_Examined", 
"Secondary_or_supplementary_reads", "Unmapped_Reads", "Unpaired_Read_Duplicates", 
"Read_Pair_Duplicates", "Read_Pair_Optical_Duplicates", "Percent_Duplicates", 
"Estimated_Library_Size")

write.table(dup.logs, "../Duplication_Metrics2.txt", sep="\t", quote=FALSE)
