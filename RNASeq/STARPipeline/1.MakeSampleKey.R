#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

# Make a sample-fastq key for use downstream

sample.key <- read.delim(args[1], stringsAsFactors=FALSE)
# This is pasted from the emails from Core for each lane of data received and looks
# like this:
# FileName	Sample Name
# WTCHG_384296_201191	S97
# WTCHG_384296_202179	S328
# WTCHG_384296_203167	S81

mapping.key <- matrix(nrow=length(unique(sample.key[, 2])), ncol=2)
mapping.key[, 1] <- unique(sample.key[, 2])
for(i in 1:nrow(mapping.key)){
  files <- c(sample.key[, 1][sample.key[, 2] == mapping.key[i, 1]])
  mapping.key[i, 2] <- paste(files, collapse=",")
}

mapping.key <- data.frame(mapping.key)

write.table(mapping.key, "mapping.info.txt", sep="\t", quote=FALSE, 
            row.names=FALSE, col.names=FALSE)
