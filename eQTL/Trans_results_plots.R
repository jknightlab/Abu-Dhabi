# trans
load("Peak_eQTL.RData")
trans.eqtl <- read.delim("Significant_trans_eQTL_Bonferroni.txt", stringsAsFactors=F)
trans.eqtl <- trans.eqtl[order(trans.eqtl$eQTL_pval), ]
trans.eqtl <- trans.eqtl[!duplicated(trans.eqtl$Gene), ]
trans.geno <- readRDS("trans.geno.rds")

eqtl.plot.side <- function(gene, snp){
  
  eqtl <- data.frame("Expression"=exprs[gene, ], 
                     "Genotype"=as.factor(trans.geno[, snp]))
  eqtl <- eqtl[complete.cases(eqtl), ]
  
  ggplot(eqtl, aes(Genotype, Expression, colour=Genotype)) +
    geom_boxplot() +
    geom_point(position=position_jitter(width=0.2)) +
    xlab(label=snp) +
    ylab(label=paste("Adjusted", gene, "Expression")) +
    theme_bw() +
    ggtitle("eQTL (Linear Model)")
}


pdf("../Trans_eQTL.pdf")
for(i in 1:nrow(trans.eqtl)){
  p <- eqtl.plot.side(trans.eqtl$Gene[i], trans.eqtl$SNP[i])
  print(p)
}
dev.off()
