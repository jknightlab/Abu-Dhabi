# To obtain transcriptional start site (TSS) information and overlay this with eQTL data

library(ggplot2)
library(data.table)
  
genelist <- read.delim("U:/Abu-Dhabi/RNASeq/Full_transcript_metadata.txt")
eqtl <- load("U:/Abu-Dhabi/eQTL/ShinyApp/Peak_eQTL.RData")
pos <- fread("U:/Abu-Dhabi/eQTL/AD_736_multi_ethnic_chip_eQTL_genotyping_b38.bim")

peak.eQTL$SNP_pos <- pos$V4[match(peak.eQTL$SNP, pos$V2)]
peak.eQTL$Gene_pos <- genelist$TSS[match(peak.eQTL$Gene, genelist$gene_name)]

peak.eQTL$Distance <- peak.eQTL$SNP_pos - peak.eQTL$Gene_pos

# Density plot distance from TSS
ggplot(peak.eQTL, aes(x=Distance)) + 
  geom_density(alpha=.3) + theme_bw() + xlab("Distance to TSS") + 
  theme(legend.position="none") + 
  geom_vline(xintercept=-250000) +
  geom_vline(xintercept=250000)

# Distance from TSS against significance
ggplot(peak.eQTL, aes(x=Distance, y=-log(eQTL_pval))) + 
  theme_bw() + 
    geom_point() + 
  geom_vline(xintercept=-250000) +
  geom_vline(xintercept=250000)

ggplot(peak.eQTL, aes(x=Distance, y=abs(eQTL_beta), colour=-log(eQTL_pval))) + 
  theme_bw() + 
  geom_point() + 
  geom_vline(xintercept=-250000) +
  geom_vline(xintercept=250000) +
  scale_color_continuous(low="darkblue", high="yellow")

ggplot(peak.eQTL, aes(x=Distance, y=eQTL_t)) + 
  theme_bw() + 
  geom_point() + 
  geom_vline(xintercept=-250000) +
  geom_vline(xintercept=250000)
