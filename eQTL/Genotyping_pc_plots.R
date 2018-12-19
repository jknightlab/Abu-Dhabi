# Genotyping PCs
library(ggplot2)

geno.pc <- read.table("eQTL/AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed_noLD_noLD_genotyping_pca_clean.eigenvec", 
                      sep="", header=F, row.names=2)
eigenvals <- read.delim("eQTL/AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed_noLD_noLD_genotyping_pca_clean.eigenval", 
                        sep="", header=F)
eigenvals$V1 <- round(eigenvals$V1, 2)

geno.pc$V1 <- NULL

pop.outliers <- read.delim("../Genotyping/pop_outliers.txt", header=F)
geno.pc$Outlier <- as.factor(rownames(geno.pc) %in% pop.outliers$V1)

for(x in seq(1, 19, 2)){
  p <- ggplot(geno.pc, aes_string(colnames(geno.pc)[x], colnames(geno.pc)[x + 1])) + 
    geom_point(aes(colour=Outlier)) + 
    theme_bw() + 
    xlab(paste("PC", x, "-", eigenvals[x, 1], "%")) + 
    ylab(paste("PC", x + 1, "-", eigenvals[x + 1, 1], "%"))  
  print(p)
  readline(prompt="Press [enter] to continue")
}

to.rm <- geno.pc[(geno.pc$V6 < -0.2 |
                   geno.pc$V6 > 0.1 |
                   geno.pc$V7 > 0.05 |
                   geno.pc$V7 < -0.12 |
                   geno.pc$V8 > 0.25 |
                   geno.pc$V9 < -0.05 |
                   geno.pc$V9 > 0.06 |
                   geno.pc$V11 < -0.25 |
                   geno.pc$V11 > 0.1 |
                   geno.pc$V12 > 0.15 |
                   geno.pc$V13 > 0.06 |
                   geno.pc$V13 < -0.07 |
                   geno.pc$V14 > 0.1 |
                   geno.pc$V14 < -0.1 |
                   geno.pc$V15 > 0.1 |
                   geno.pc$V15 < -0.8 |
                   geno.pc$V16 < -0.1 |
                   geno.pc$V16 > 0.1 |
                   geno.pc$V17 > 0.1 |
                   geno.pc$V17 < -0.1 |
                   geno.pc$V18 < -0.12 |
                   geno.pc$V18 > 0.1 |
                   geno.pc$V19 > 0.2 |
                   geno.pc$V19 < -0.1 |
                   geno.pc$V20 < -0.1 |
                   geno.pc$V20 > 0.1 |
                   geno.pc$V21 > 0.2 |
                   geno.pc$V22 < -0.12 |
                   geno.pc$V22 > 0.1), ]

geno.pc$Outlier <- as.factor(rownames(geno.pc) %in% rownames(to.rm))

pdf("Genotyping_PCs_pop_outliers.pdf", onefile=T, useDingbats=F)
for(x in seq(1, 19, 2)){
  p <- ggplot(geno.pc, aes_string(colnames(geno.pc)[x], colnames(geno.pc)[x + 1])) + 
    geom_point(aes(colour=Outlier)) + 
    theme_bw() + 
    xlab(paste("PC", x, "-", eigenvals[x, 1], "%")) + 
    ylab(paste("PC", x + 1, "-", eigenvals[x + 1, 1], "%"))  
  print(p)
  readline(prompt="Press [enter] to continue")
}
dev.off()

ggplot(eigenvals, aes(1:20, V1)) + geom_col() + theme_bw() + geom_hline(yintercept=1)

pop.outliers <- rownames(to.rm)
write.table(pop.outliers, "pop_outliers_Katie.txt", sep="\t", quote=F, 
            row.names=F, col.names=F)








geno.pc <- read.table("../Genotyping/AD_736_multi_ethnic_chip_updated_inds_snps_removed_noLD_noLD_genotyping_pca_clean.eigenvec", sep="", header=F, row.names=2)
eigenvals <- read.delim("../Genotyping/AD_736_multi_ethnic_chip_updated_inds_snps_removed_noLD_noLD_genotyping_pca_clean.eigenval", 
                        sep="", header=F)
eigenvals$V1 <- round(eigenvals$V1, 2)

geno.pc$V1 <- NULL

for(x in seq(1, 19, 2)){
  p <- ggplot(geno.pc, aes_string(colnames(geno.pc)[x], colnames(geno.pc)[x + 1])) + 
    geom_point() + 
    theme_bw() + 
    xlab(paste("PC", x, "-", eigenvals[x, 1], "%")) + 
    ylab(paste("PC", x + 1, "-", eigenvals[x + 1, 1], "%"))  
  print(p)
  readline(prompt="Press [enter] to continue")
}
