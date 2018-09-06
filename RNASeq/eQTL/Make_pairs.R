library(data.table)

gene.info <- read.delim("Full_transcript_metadata.txt")

genes <- gene.info[, c("chr", "TSS")]
genes$chr <- gsub("chr", "", genes$chr)

dim(gene.info[which(is.na(gene.info$TSS)),])
#There are 28 genes without TSS information 

genes.list <- list()
for(i in 1:22){
  genes.list[[i]] <- subset(genes, chr == i)
}

# dataframe with SNP names as row names and columns for chromosome number and SNP position
# for each chr
bim <- data.frame(fread(paste("eQTL/AD_736_multi_ethnic_chip_update_strand_inds_snps_removed_no_parents.bim", 
                              sep=""), header=FALSE))

chr.bim.list <- list()
for(i in 1:22){
  chr.bim.list[[i]] <- subset(bim, V1 == i)
  rownames(chr.bim.list[[i]]) <- as.character(chr.bim.list[[i]][, 2])
  chr.bim.list[[i]] <- chr.bim.list[[i]][, c(1, 4)]
}

for(c in 1:22){
  results=NULL
  for (i in 1:nrow(genes.list[[c]])) {
    snp.names <- rownames(chr.bim.list[[c]])[which(as.numeric(genes.list[[c]][i, 1]) == chr.bim.list[[c]][, 1] 
                                                   & (chr.bim.list[[c]][, 2] > (genes.list[[c]][i, 2] - 250000) 
                                                      &  chr.bim.list[[c]][, 2] < (genes.list[[c]][i, 2] + 250000)))]
    if(length(snp.names) > 1){
      gene.snp.pair <- cbind(rownames(genes.list[[c]])[i], snp.names)
      results <- rbind(results, gene.snp.pair)
    }
    
  }
  
  results.2 <- data.frame(results)
  write.table(results.2, paste("eQTL/chr", c, "Gene_snp_pairs.txt", sep=""), sep="\t", 
              quote=FALSE, row.names=FALSE)
}
