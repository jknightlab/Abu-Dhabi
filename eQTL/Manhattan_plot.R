# eQTL plots

library(data.table)
library(qqman)
library(RColorBrewer)

eqtl.results <- base::readRDS("/well/jknight/AbuDhabiRNA/eQTL/Results/cis/eQTL_results_final_lm_cis.rds")
sig.results <- read.delim("/well/jknight/AbuDhabiRNA/eQTL/Results/cis/All_Sig_Results_eQTL.txt")

# Manhattan plot
# on galahad
bim <- data.frame(fread("/well/jknight/AbuDhabiRNA/eQTL/Genotyping/AD_736_multi_ethnic_chip_clean_b38.bim",
                                         header=FALSE))
snp.pos <- bim[, c(2, 1, 4)]
colnames(snp.pos) <- c("SNP", "Chrm_snp", "Pos")
snp.pos$SNP <- gsub(":", ".", snp.pos$SNP)

results.for.manhattan <- cbind(SNP=as.character(eqtl.results$SNP),
                               CHR=snp.pos$Chrm_snp[match(eqtl.results$SNP,
                                                          snp.pos$SNP)],
                               BP=snp.pos$Pos[match(eqtl.results$SNP,
                                                    snp.pos$SNP)],
                               P=eqtl.results$eQTL_pval)

results.for.manhattan <- data.frame(results.for.manhattan)
results.for.manhattan$CHR <- as.integer(as.character(results.for.manhattan$CHR))
results.for.manhattan$BP <- as.integer(as.character(results.for.manhattan$BP))
results.for.manhattan$P <- as.numeric(as.character(results.for.manhattan$P))

results.for.manhattan$Sig <- paste0(as.character(eqtl.results$SNP), as.character(eqtl.results$Gene)) %in%
  paste0(as.character(sig.results$SNP), as.character(sig.results$Gene))

results.for.manhattan <- results.for.manhattan[complete.cases(results.for.manhattan), ]

png("eQTLManhattan.png", width=960)
manhattan(results.for.manhattan, genomewideline=FALSE,
          suggestiveline=FALSE, 
          col=c(brewer.pal(9, "Set1"), brewer.pal(12, "Set3"),
                brewer.pal(8, "Set2")))
dev.off()

png("eQTLManhattanSig.png", width=960)
manhattan(results.for.manhattan, genomewideline=FALSE,
          suggestiveline=FALSE, highlight=results.for.manhattan$SNP[results.for.manhattan$Sig == TRUE])
dev.off()