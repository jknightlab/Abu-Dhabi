#!/usr/bin/env Rscript

# Mapping QC
library(ggplot2)

log.paths <- list.files(path=".", pattern="*Log.final.out", full.names=TRUE, recursive=TRUE)
logs <- lapply(log.paths, read.delim, row.names=1, header=FALSE)
logs <- lapply(logs, as.data.frame)

cbindlist <- function(list) {
  n <- length(list)
  res <- list[[1]]
  for (i in 2:n) res <- cbind(res, list[[i]])
  return(res)
}

mapping.stats <- cbindlist(logs)
mapping.stats <- t(mapping.stats)
colnames(mapping.stats) <- gsub(" ", "", colnames(mapping.stats))
colnames(mapping.stats) <- gsub("\\|", "", colnames(mapping.stats))
mapping.stats <- data.frame(mapping.stats)

pdf("MappingQC.pdf", onefile=TRUE)

ggplot(mapping.stats, aes(Numberofinputreads, Uniquelymappedreadsnumber)) +
  geom_point() +
  theme_bw() +
  xlab(label="Number of reads") +
  ylab(label="Number of uniquely mapping reads")

ggplot(mapping.stats, aes(Numberofinputreads, Numberofsplices.Total)) +
  geom_point() +
  theme_bw() +
  xlab(label="Number of reads") +
  ylab(label="Number of splices identified")

ggplot(mapping.stats, aes(Numberofinputreads, Mismatchrateperbase..)) +
  geom_point() +
  theme_bw() +
  xlab(label="Number of reads") +
  ylab(label="Mismatch rate per base")

ggplot(mapping.stats, aes(Numberofinputreads, Numberofreadsmappedtomultipleloci)) +
  geom_point() +
  theme_bw() +
  xlab(label="Number of reads") +
  ylab(label="Number of multi-mapping reads")

dev.off()
