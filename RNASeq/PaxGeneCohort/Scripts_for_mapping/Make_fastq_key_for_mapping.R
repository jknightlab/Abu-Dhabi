# Make sample-fastq key for mapping

library(gdata)

sample.key <- read.delim("~/Abu-Dhabi/RNASeq/PaxGeneCohort/AD_full_sample_fastq_key.txt",
                         stringsAsFactors=FALSE)
core.summary <- read.xls("~/Abu-Dhabi/RNASeq/PaxGeneCohort/P170245_report_15Mar18.xlsx", 
                         sheet=3, stringsAsFactors=FALSE)

table(core.summary$Sample.Name %in% sample.key$Sample.Name, core.summary$Sequenced.library.)
# No Yes
# FALSE   5   0
# TRUE   13 623
# fastq files for the 623 successful samples, and 5 that failed

# Make a list of the 623 samples and the corresponding fastq files
# NB FG84M is duplicated so was deleted manually
mapping.key <- matrix(nrow=623, ncol=2)
mapping.key[, 1] <- core.summary$Sample.Name[core.summary$Sequenced.library. == "Yes"]
for(i in 1:nrow(mapping.key)){
  files <- c(sample.key$FileName[sample.key$Sample.Name == mapping.key[i, 1]])
  mapping.key[i, 2] <- paste(files, collapse=",")
}

mapping.key <- data.frame(mapping.key)

write.table(mapping.key, "/well/jknight/AbuDhabiRNA/Katie/mapping/2.mapping/mapping.info.txt",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)



################################################################################

# Make index of bam files for mapping QC

s.info <- read.delim("/well/jknight/AbuDhabiRNA/Katie/mapping/2.mapping/mapping.info.txt")

bam.files <- data.frame("Sample ID"=s.info[, 1], 
                        "Bam File"=paste("/well/jknight/AbuDhabiRNA/Katie/mapping/2.mapping/", 
                                         s.info[, 1], "/accepted_hits.bam", sep=""), 
                        "Notes"=NA)
write.table(bam.files, "QC/bam.files.txt", sep="\t", quote=FALSE)
