library(data.table)

gene.info <- read.delim("Full_transcript_metadata.txt")

genes <- gene.info[, c("chr", "TSS")]
rownames(genes) <- gene.info$gene_name
genes$chr <- gsub("chr", "", genes$chr)

dim(gene.info[which(is.na(gene.info$TSS)),])

bim <- data.frame(fread(paste("eQTL/AD_736_multi_ethnic_chip_eQTL_genotyping_b38.bim", 
                              sep=""), header=FALSE))
rownames(bim) <- as.character(bim[, 2])
bim <- bim[, c(1, 4)]

results=NULL

irange <- 1:nrow(genes)
results <- do.call(rbind, lapply(irange, function(i){
  snp.names <- rownames(bim)[which(as.numeric(genes[i, 1]) == bim[, 1]
                                   & (bim[, 2] > (genes[i, 2] - 250000) 
                                      &  bim[, 2] < (genes[i, 2] + 250000)))]
  if(length(snp.names) > 1){
    gene.snp.pair <- cbind(rownames(genes)[i], snp.names)
    
  }
    
}))
  
results.2 <- data.frame(results)
write.table(results.2, "eQTL/Gene_snp_pairs.txt", sep="\t", quote=FALSE, row.names=FALSE)
