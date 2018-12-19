library(data.table)

exprs <- read.delim("~/Abu-Dhabi/RNASeq/eQTL/Normalised_GEX_data.txt")
colnames(exprs) <- gsub("X", "", colnames(exprs))
exprs <- as.matrix(exprs)

info <- read.delim("~/Abu-Dhabi/RNASeq/eQTL/sample_info.txt")
rownames(info) <- info$SampleID

gex.pcs <- read.delim("~/Abu-Dhabi/RNASeq/eQTL/GEX_PCs.txt")
gex.pcs <- as.matrix(gex.pcs)

geno.pcs <- read.delim("~/Abu-Dhabi/Genotyping/AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed_noLD_noLD_genotyping_pca_clean.eigenvec", sep="", row.names=2, header=F)
geno.pcs$V1 <- NULL
geno.pcs <- as.matrix(geno.pcs)

geno <- data.frame(fread("/well/jknight/AbuDhabiRNA/eQTL/Genotyping/AD_736_multi_ethnic_chip_eQTL_genotyping_b38.raw"))
rownames(geno) <- geno[, 1]
geno[, 1:6] <- NULL
colnames(geno) <- gsub("X", "", colnames(geno))
colnames(geno) <- substr(colnames(geno), 1, nchar(colnames(geno))-2)
geno <- as.matrix(geno)
  
pairs <- read.delim("~/Abu-Dhabi/RNASeq/eQTL/Gene_snp_pairs.txt")
pairs[, 1] <- as.character(pairs[, 1])
pairs[, 2] <- make.names(pairs[, 2])
pairs[, 2] <- gsub("X", "", pairs[, 2])
colnames(pairs) <- c("Gene", "SNP")
  
PCs <- cbind(gex.pcs[, 1:25], geno.pcs[, 1:4])
expression.set <- exprs
expression.set <- t(expression.set)
expression.pc <- PCs
num.PC <- ncol(expression.pc)
  
# regress out expression
regressed.data <- matrix(nrow=nrow(expression.set), ncol=ncol(expression.set))
for(i in 1:(dim(expression.set)[2])){
  model <- lm(as.matrix(expression.set[, i]) ~ as.matrix(expression.pc[, 1:(num.PC)]))
  regressed.data[, i] <- expression.set[, i] - 
    rowSums(sapply(1:(num.PC), function(i)model$coefficients[i+1]*expression.pc[, 1:(num.PC)][,i]))
}
rownames(regressed.data) <- rownames(expression.set)
colnames(regressed.data) <- colnames(expression.set)
regressed.data <- t(regressed.data)

exprs <- regressed.data

save(list=c("exprs", "geno", "info", "gex.pcs", "geno.pcs", "pairs"),
       file = "/well/jknight/AbuDhabiRNA/eQTL/eQTL.25PCs.RData")
  
