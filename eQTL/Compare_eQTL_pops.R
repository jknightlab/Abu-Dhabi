# Compare eQTL results

library(ggplot2)
library(data.table)


# HV
hv <- read.delim("U:/eQTL/NESDA_NTR_Conditional_eQTL_Catalog/eqtl_base/cis_eQTL_table_conditional_ALL",
                 stringsAsFactors=F)
hv$Alt <- substr(hv$SNP, nchar(hv$SNP), nchar(hv$SNP))
hv <- hv[hv$Alt %in% c("A", "C", "T", "G"), ]
hv$GeneSNP <- paste0(hv$Gene, hv$rs.SNP)
genes <- read.delim("U:/eQTL/NESDA_NTR_Conditional_eQTL_Catalog/genes_tested/genes_tested.txt",
                    header=F)
snps <- read.delim("U:/eQTL/NESDA_NTR_Conditional_eQTL_Catalog/SNPs_tested_RSID/SNPs_tested_RSID",
                   header=F)

# T2D
t2d <- readRDS("U:/Abu-Dhabi/eQTL/ShinyApp/eQTL_results_final_lm_cis.rds")
t2d <- t2d[t2d$Gene %in% genes$V1, ]
t2d <- t2d[t2d$SNP %in% snps$V1, ]
t2d$GeneSNP <- paste0(t2d$Gene, t2d$SNP)
bim <- fread("../eQTL/AD_736_multi_ethnic_chip_eQTL_genotyping_b38.bim")
t2d$Alt <- bim$V5[match(t2d$SNP, bim$V2)]
t2d$Ref <- bim$V6[match(t2d$SNP, bim$V2)]

t2d.sig <- read.delim("U:/Abu-Dhabi/eQTL/All_Sig_Results_eQTL.txt", stringsAsFactors=F)


m <- merge(hv, t2d, by="GeneSNP", all.y=T)
m$Beta[is.na(m$Beta)] <- 0

ggplot(m, aes(Beta, eQTL_beta)) + 
  geom_point() +
  theme_bw()



m$action <- "NA"
m[which(m$Alt.x == m$Alt.y), "action"] <- "same"

table(m$action)
# NA  same 
# 71624 59774 

# Alleles can't distinguish because of strand
m[which(m$Alt.y == "A" & m$Ref == "T" |
                     m$Alt.y == "T" & m$Ref == "A" |
                     m$Alt.y == "C" & m$Ref == "G" |
                     m$Alt.y == "G" & m$Ref == "C"), 
             "action"] <-"strand"

table(m$action)
# NA   same strand 
# 67951  55946   7501

# Same but different strand
m[which(m$Alt.x == "A" & m$Alt.y == "T" & 
                     m$action != "strand" |
                     m$Alt.x=="T" & m$Alt.y=="A" & 
                     m$action != "strand" |
                     m$Alt.x=="C" & m$Alt.y=="G" & 
                     m$action != "strand" |
                     m$Alt.x=="G" & m$Alt.y=="C" & 
                     m$action != "strand"),
             "action"] <-"same"

table(m$action)
# NA   same strand 
# 11341 112556   7501  

m[which(m$action == "NA"), "action"] <- "flip"

table(m$action)
# flip   same strand 
# 11341 112556   7501 

m$eQTL_beta.2 <- m$eQTL_beta
m[which(m$action == "flip"), "eQTL_beta.2"] <- 
  m[which(m$action == "flip"), "eQTL_beta.2"] * (-1)

# Remove 319 SNPs which can't be distinguished because of strand
m$eQTL_beta.2[which(m$action == "strand")] <- NA
m <- m[complete.cases(m), ]

ggplot(m, aes(Beta, eQTL_beta.2)) + 
  geom_point() +
  theme_bw()



hv.red <- hv[order(hv$P), ]
hv.red <- hv.red[!duplicated(hv.red$Gene), ]

t2d.red <- t2d.sig[order(t2d.sig$eQTL_pval), ]
t2d.red <- t2d.red[!duplicated(t2d.red$Gene), ]

m.red <- merge(hv.red, t2d.red, by="Gene")
plot(order(m.red$P), order(m.red$eQTL_pval))
plot(m.red$Beta, m.red$eQTL_beta)
