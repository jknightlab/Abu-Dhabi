# GWAS overlap

#################################################################################
#
# Colocalisation of GWAS and eQTL hits
#
#################################################################################

# AIM: to find loci that are significant in both the sepsis GWAS and eQTL
#      test for a shared causal effect

# STRATEGY:
# Start with the GWAS summary stats
# Find all that are genome-wide significant (5e-08)
# Export this list
# Get all SNPs in LD (r2>0.8) from ?1K genomes
# Look for overlap with eQTL (all significant)
# Taking the GWAS hits with an overlap
# Loop through all the GWAS outputs to define independent loci
# Find all eQTL SNPs in same region and in LD (r2>0.5)
# Make local assn plots

# To test for shared effect with coloc:
# Reduce to SNPs tested in both
# Look up eQTL results for the same SNPs
# Pull out MAF etc for the GWAS to calculate posterior probabliities
# Run coloc on each locus

options(stringsAsFactors = FALSE)

library(coloc)
library(rlist)
library(gdata)
library(survival)
library(ggplot2)

# read in results
setwd("U:/Abu-Dhabi/eQTL/")
GWAS <- read.table("../GWAS/AD2463_T2D_covar.assoc.logistic", header=TRUE)
load("ShinyApp/Peak_eQTL.RData")
eqtl.results <- all.eqtl
rm(all.eqtl)
bim <- read.delim("../GWAS/AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed.bim", header=F)
eqtl.results$CHR <- bim$V1[match(eqtl.results$SNP, bim$V2)]
eqtl.results$BP <- bim$V4[match(eqtl.results$SNP, bim$V2)]

# MAF
maf <- read.table("../GWAS/AD2463.frq", header=T)
GWAS$MAF <- maf$MAF[match(GWAS$SNP, maf$SNP)]

# Significant GWAS SNPs
GWAS <- GWAS[order(GWAS$P), ]
GWAS$eQTL <- GWAS$SNP %in% eqtl.results$SNP

GWAS$SNPFull <- paste0(GWAS$CHR, ".", GWAS$BP)
eqtl.results$SNPFull <- paste0(eqtl.results$CHR, ".", eqtl.results$BP)

# export SNPs at suggestive significance
gwas.hits <- subset(GWAS, P < 5e-5)
combine <- merge(gwas.hits, eqtl.results, by.x="SNPFull", by.y="SNPFull")
# 
# combine.all <- merge(GWAS, eqtl.results, by.x="SNP", by.y="SNP")
# combine.all <- combine.all[complete.cases(combine.all), ]
# ggplot(combine.all, aes(P, eQTL_pval)) +
#   geom_point() + theme_bw() + geom_vline(xintercept=5e-5)
# 
# # Manhattan plot
# library(dplyr)
# library(ggrepel)
# 
# snpsOfInterest <- as.character(peak.eQTL$SNP)
# don <- GWAS %>%
# 
#   filter(-log10(P)>1) %>%
#   # Compute chromosome size
#   group_by(CHR) %>%
#   summarise(chr_len=max(BP)) %>%
# 
#   # Calculate cumulative position of each chromosome
#   mutate(tot=cumsum(chr_len)-chr_len) %>%
#   select(-chr_len) %>%
# 
#   # Add this info to the initial dataset
#   left_join(GWAS, ., by=c("CHR"="CHR")) %>%
# 
#   # Add a cumulative position of each SNP
#   arrange(CHR, BP) %>%
#   mutate( BPcum=BP+tot) %>%
# 
#   # Add highlight and annotation information
#   mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
#   mutate( is_annotate=ifelse(-log10(P)>7, "yes", "no"))
# 
# axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
# 
# png("GWAS_Manhattan_peak_eQTL.png", width=1400)
# print(ggplot(don, aes(x=BPcum, y=-log10(P))) +
# 
#   # Show all points
#   geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
#   scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
# 
#   # custom X axis:
#   scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
#   scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
# 
#   # Add highlighted points
#   geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
# 
#   # Add label using ggrepel to avoid overlapping
#   #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
# 
#   # Custom the theme:
#   theme_bw() +
#   theme(
#     legend.position="none",
#     panel.border = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank()
#   ))
# dev.off()

################################################################################
# for each line in the file extract ld info
snps <- bim[bim$V2 %in% gwas.hits$SNP, ]
write.table(snps, "../GWAS/GWAS_eQTL_snps_to_extract.txt", row.names=F)
# dos2unix ~/Abu-Dhabi/GWAS/GWAS_eQTL_snps_to_extract.txt
#
# while read -r line
# do
#   chr="$(echo "$line" | cut -f 1)"
#   bp="$(echo "$line" | cut -f 4)"
#   snp="$(echo "$line" | cut -f 2)"
#
#   ~/plink --bfile /well/jknight/Sepsis_eQTL/1kg_EUR/1000genomes_EUR_chr$chr \
#   --chr $chr \
#   --from-bp $bp \
#   --to-bp $bp \
#   --make-bed \
#   --out ~/eQTL/GWAS\ overlap/1kgLD/tmp.$snp
#
#   read -r ids < ~/eQTL/GWAS\ overlap/1kgLD/tmp.$snp.bim
#   rsid="$(echo "$ids" | cut -f 2)"
#   rm ~/eQTL/GWAS\ overlap/1kgLD/tmp*
#
#   ~/plink --bfile /well/jknight/Sepsis_eQTL/1kg_EUR/1000genomes_EUR_chr$chr \
#   --ld-snp $rsid \
#   --r2 \
#   --ld-window 99999 \
#   --ld-window-kb 3000 \
#   --ld-window-r2 0 \
#   --out ~/eQTL/GWAS\ overlap/1kgLD/$snp
#
# done < ~/eQTL/GWAS\ overlap/GenOSept_suggestive_snps_to_extract.txt

# rm ~/eQTL/GWAS\ overlap/1kgLD/*.log
# rm ~/eQTL/GWAS\ overlap/1kgLD/*.nosex
#################################################################################

# LD.files <- list.files(path="U:/eQTL/GWAS overlap/1kgLD/", "*.ld",
#                        full.names=TRUE)
# LD.stats <- lapply(LD.files, read.table, header=TRUE)
# names(LD.stats) <- list.files(path="U:/eQTL/GWAS overlap/1kgLD/", "*.ld",
#                               full.names=FALSE)
#
# # keep snps in ld with gwas snp
# for(i in 1:length(LD.stats)){
#   snp <- names(LD.stats)[i]
#   snp <- substr(snp, 1, nchar(snp)-7)
#   snp <- as.numeric(unlist(strsplit(snp, split="[.]")))
#   snp <- which(LD.stats[[i]]$CHR_A == snp[1] &
#                  LD.stats[[i]]$BP_A == snp[2] &
#                  LD.stats[[i]]$R2 >0.8)
#   LD.stats[[i]] <- LD.stats[[i]][snp, ]
# }
#
# # pull out these snps from the gwas and eqtl results
# combine <- list()
# for(i in 1:length(LD.stats)){
#   tryCatch({
#     combine[[i]] <- LD.stats[[i]][, 4:7]
#     colnames(combine[[i]]) <- c("CHR", "BP", "SNP", "R2")
#     combine[[i]]$SNP.sh <- paste0(combine[[i]]$CHR, ".", combine[[i]]$BP)
#     combine[[i]] <- merge(combine[[i]], sig.eqtl, by.x="SNP.sh", by.y="SNP.sh")
#   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
#
# }
# n1 <- unlist(lapply(combine, nrow))
#
# # One locus
# combine <- combine[[which(n1 > 0)]]
combine

# pull out all gwas snps in window and LD
merging.distance <- 1000000
r2.threshold <- 0.5

gwas.locus <- data.frame(GenOSept[(GenOSept$CHR == combine$CHR[1] &
                                     GenOSept$BP %in% (combine$BP[1] - merging.distance):(combine$BP[1] + merging.distance)), ])

gwas.locus <- merge(gwas.locus, eqtl.results[eqtl.results$Symbol == "NCOA4", ],
                    by.x="SNP", by.y="SNP.sh")

ld <- read.table("../GWAS overlap/rs11004438.ld", header=TRUE)
ld <- subset(ld, R2 >=0.5)
gwas.locus <- gwas.locus[gwas.locus$rsID %in% ld$SNP_B |
                           gwas.locus$rsID %in% ld$SNP_A, ]

gwas.maf <- read.delim("../GWAS overlap/GAinS_chr10_locus_MAF.txt")
gwas.maf$rsid <- gsub("-", ".", gwas.maf$rsid)
gwas.locus$MAF <- gwas.maf$all_maf[match(gwas.locus$SNP, gwas.maf$rsid)]

gwas <- list(pvalues=gwas.locus$P,
             N=2534,
             MAF=gwas.locus$MAF,
             type="cc",
             s=0.226)
eqtl <- list(pvalues=gwas.locus$eQTL_pval,
             N=592,
             MAF=snp.name.key$MAF[match(gwas.locus$SNP, snp.name.key$SNP)],
             beta=gwas.locus$eQTL_beta,
             varbeta=gwas.locus$eQTL_SE^2,
             type="quant", sdY=0.3937727)

results <- coloc.abf(gwas, eqtl)

# PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf
# 7.63e-24  3.41e-24  1.15e-01  5.06e-02  8.34e-01
# [1] "PP abf for shared variant: 83.4%"

library(kburnhamfunctions)

# GWAS local assn plot
ld <- read.table("../GWAS overlap/rs11004438.ld", header=TRUE)
gwas.locus <- data.frame(GenOSept[(GenOSept$CHR == combine$CHR[1] &
                                     GenOSept$BP %in% (combine$BP[1] - merging.distance):(combine$BP[1] + merging.distance)), ])
gwas.locus$rsID <- snp.name.key$rsID[match(gwas.locus$SNP, snp.name.key$SNP)]

gwas.result <- data.frame(SNP=gwas.locus$rsID,
                          POS=gwas.locus$BP,
                          PVAL=gwas.locus$P,
                          TYPE="typed")
gwas.result <- gwas.result[complete.cases(gwas.result), ]
rownames(gwas.result) <- gwas.result$SNP
gwas.result$RSQR <- ld$R2[match(gwas.result$SNP, ld$SNP_B)]

omni <- read.delim("U:/GAinS/Genotyping/GAinS_genotyping_275samples.bim",
                   header=FALSE)
exome <- read.delim("U:/GAinS/Genotyping/clean_exome_illuminus.bim",
                    header=FALSE)
# gwas.result$TYPE[!(gwas.result$SNP %in% omni$V2 &
# gwas.result$SNP %in% exome$V2)] <- "imputed"
gwas.result <- gwas.result[complete.cases(gwas.result), ]
make.fancy.locus.plot("rs11004435", "NCOA4 region", "10", gwas.result, 5, 4.314e-05)

# eQTL local assn plot
# Pull out all results for same region
eqtl.results$CHR <- snp.name.key$Chr[match(eqtl.results$SNP.sh, snp.name.key$SNP)]
eqtl.results$BP <- snp.name.key$BP[match(eqtl.results$SNP.sh, snp.name.key$SNP)]
eqtl.result <- eqtl.results[eqtl.results$CHR == 10 &
                              eqtl.results$Symbol == "NCOA4" &
                              (eqtl.results$BP %in% (combine$BP[1] - merging.distance):(combine$BP[1] + merging.distance)), ]
eqtl.result <- data.frame(SNP=eqtl.result$rsID,
                          POS=eqtl.result$BP,
                          PVAL=eqtl.result$eQTL_pval,
                          TYPE="typed")
eqtl.result <- eqtl.result[complete.cases(eqtl.result), ]
rownames(eqtl.result) <- eqtl.result$SNP
eqtl.result$RSQR <- ld$R2[match(eqtl.result$SNP, ld$SNP_B)]
eqtl.result$TYPE[!(eqtl.result$SNP %in% omni$V2 &
                     eqtl.result$SNP %in% exome$V2)] <- "imputed"

eqtl.result <- eqtl.result[complete.cases(eqtl.result), ]
make.fancy.locus.plot("rs2125771", "NCOA4 region", "10", eqtl.result, 25, 1.423745e-25)

gwas.result <- gwas.result[gwas.result$SNP %in% eqtl.result$SNP, ]
make.fancy.locus.plot("rs11004435", "NCOA4 region", "10", gwas.result, 5, 4.314e-05)

eqtl.result <- eqtl.result[eqtl.result$POS < 51700000, ]
make.fancy.locus.plot("rs2125771", "NCOA4 region", "10", eqtl.result, 25, 1.423745e-25)
gwas.result <- gwas.result[gwas.result$SNP %in% eqtl.result$SNP, ]
make.fancy.locus.plot("rs11004435", "NCOA4 region", "10", gwas.result, 5, 4.314e-05)

# Replication of GWAS association
gains.geno <- read.table("../GWAS overlap/rs11004435_genotypes_GAinS.raw", header=T)
genosept.samples <- read.delim("N:/jknight/Data/GAinS/GenOSept_Database/Files for Access Import/All_samples_genotyped_in_GWAS.txt", header=F)

dups <- which(gains.geno$FID %in% genosept.samples$V1)
first.cohort <- 1:275
first.cohort <- first.cohort[!first.cohort %in% dups]
second.cohort <- 276:910
second.cohort <- second.cohort[!second.cohort %in% dups]

gains.geno <- gains.geno[!(gains.geno$FID %in% genosept.samples$V1), ]
survival <- read.xls("N:/jknight/Data/GAinS/Clinical data/eCRF June 2016/OUT_02jun2016.xls")
survival$SubjectBarCode <- gsub("uk", "UK", survival$SubjectBarCode)
gains.geno$FID <- gsub("uk", "UK", gains.geno$FID)
gains.geno <- gains.geno[gains.geno$FID %in% survival$SubjectBarCode, ]
time <- as.numeric(as.character(survival$DaystoDeath[match(gains.geno$FID,
                                                           survival$SubjectBarCode)]))
diagnosis <- survival$diagnosis[match(gains.geno$FID, survival$SubjectBarCode)]
time[which(time > 28)] <- NA
status <- as.numeric(!is.na(time))
time[which(is.na(time))] <- 28
genotype <- gains.geno$X10.51498493_C_A_A

# Plotting survival curves
my.survfit <- survfit(Surv(time, status) ~ genotype)

plot(my.survfit, col=c(1:3), xlab="Survival in days",
     ylab="Proportion surviving", conf.int=FALSE)
legend("bottomright", c("rs11004435 AA", "rs11004435 AC", "rs11004435 CC"),
       lty=c(1, 1), col=c(1, 2, 3))

survdiff(Surv(time, status) ~ genotype)
# p=0.05

# As 2 replication cohorts:
my.survfit <- survfit(Surv(time[first.cohort], status[first.cohort]) ~ genotype[first.cohort])
plot(my.survfit, col=c(1:3), xlab="Survival in days",
     ylab="Proportion surviving", conf.int=FALSE)
legend("bottomright", c("rs11004435 AA", "rs11004435 AC", "rs11004435 CC"),
       lty=c(1, 1), col=c(1, 2, 3))
survdiff(Surv(time[first.cohort], status[first.cohort]) ~ genotype[first.cohort])
# p=0.008

my.survfit <- survfit(Surv(time[second.cohort], status[second.cohort]) ~ genotype[second.cohort])
plot(my.survfit, col=c(1:3), xlab="Survival in days",
     ylab="Proportion surviving", conf.int=FALSE)
legend("bottomright", c("rs11004435 AA", "rs11004435 AC", "rs11004435 CC"),
       lty=c(1, 1), col=c(1, 2, 3))
survdiff(Surv(time[second.cohort], status[second.cohort]) ~ genotype[second.cohort])
# p=0.5

my.survfit <- survfit(Surv(time, status) ~ diagnosis + genotype)
plot(my.survfit, col=c(1:6), xlab="Survival in days",
     ylab="Proportion surviving", conf.int=FALSE)
legend("bottomright", c("CAP AA", "CAP AC", "CAP CC", "FP AA", "FP AC", "FP CC"),
       lty=1, col=c(1:6))
survdiff(Surv(time, status) ~ diagnosis + genotype)
# p=0.004

# Sex effect?
demo <- read.xls("N:/jknight/Data/GAinS/Clinical data/eCRF June 2016/DEMO_02jun2016.xls")
sex <- demo$sex[match(gains.geno$FID, demo$SubjectBarCode)]

survdiff(Surv(time, status) ~ genotype + sex)
# p=0.2

male <- sex == 1
my.survfit <- survfit(Surv(time[male], status[male]) ~ genotype[male])

plot(my.survfit, col=c(1:3), xlab="Survival in days (Male only)",
     ylab="Proportion surviving", conf.int=FALSE)
legend("bottomright", c("rs11004435 AA", "rs11004435 AC", "rs11004435 CC"),
       lty=c(1, 1), col=c(1, 2, 3))

survdiff(Surv(time[male], status[male]) ~ genotype[male])
# p=011

female <- sex == 2
my.survfit <- survfit(Surv(time[female], status[female]) ~ genotype[female])

plot(my.survfit, col=c(1:3), xlab="Survival in days (female only)",
     ylab="Proportion surviving", conf.int=FALSE)
legend("bottomright", c("rs11004435 AA", "rs11004435 AC", "rs11004435 CC"),
       lty=c(1, 1), col=c(1, 2, 3))

survdiff(Surv(time[female], status[female]) ~ genotype[female])
# p= 0.3

#################################################################################

all.exprs <- read.delim("C:/Users/kburnham/Desktop/Raw_Gains_Data/All_GEX_dedup_filt_norm_combat.txt")
# pcs <- prcomp(t(all.exprs))$x
# resids <- matrix(NA, ncol = ncol(all.exprs), nrow = nrow(all.exprs))
# rownames(resids) <- rownames(all.exprs)
# colnames(resids) <- colnames(all.exprs)
#
# for(i in 1:nrow(all.exprs)){
#   data <- as.data.frame(cbind(t(all.exprs[i, ]), pcs[, 1:30]))
#   colnames(data) <- c('Expression', colnames(pcs)[1:30])
#   model <- lm(Expression ~ ., data = data)
#   resids[i, ] <- model$residuals
# }
# all.exprs <- data.frame(resids)

all.info <- read.table("C:/Users/kburnham/Desktop/Raw_Gains_Data/SRS_info.txt")
all.info$Survival <- is.na(all.info$DaystoDeath)
all.info$Survival[as.numeric(all.info$DaystoDeath) > 28] <- TRUE

first.info <- all.info[order(all.info$Day, decreasing=T), ]
first.info <- first.info[!duplicated(first.info$GAinSID), ]
first.exprs <- all.exprs[, match(first.info$SampleID, colnames(all.exprs))]

de.boxplot(first.exprs[which(rownames(first.exprs) == "1300671"), ],
           first.info$Survival, cov.name="Survival", gene.name="NCOA4")

de.boxplot(first.exprs[which(rownames(first.exprs) == "1300671"), ],
           first.info$Sex, cov.name="Sex", gene.name="NCOA4")

de.boxplot(first.exprs[which(rownames(first.exprs) == "1300671"), ],
           first.info$Diagnosis, cov.name="Diagnosis", gene.name="NCOA4")

de.boxplot(all.exprs[which(rownames(all.exprs) == "1300671"),
                     all.info$SRS_New == 1],
           all.info$Survival[all.info$SRS_New == 1],
           cov.name="Survival", gene.name="NCOA4")

de.boxplot(all.exprs[which(rownames(all.exprs) == "1300671"),
                     all.info$SRS_New == 2],
           all.info$Survival[all.info$SRS_New == 2],
           cov.name="Survival", gene.name="NCOA4")

de.boxplot(all.exprs[which(rownames(all.exprs) == "1300671"), ],
           covariate=all.info$SRS_New,
           cov.name="SRS", gene.name="NCOA4")

de.boxplot(first.exprs[which(rownames(first.exprs) == "1300671"), ],
           covariate=first.info$SRS_New, col.group=first.info$Survival,
           cov.name="SRS", gene.name="NCOA4")










































# pull out all overlapping eQTL results

# GenOSept$N.eQTL <- NA
# for(i in 1:nrow(GenOSept)){
#   n <- length(which(sig.eqtl$SNP == GenOSept$SNP_full[i]))
#   GenOSept$N.eQTL[i] <- n
# }

# gwas.loci <- data.frame(gwas.hits[1, ])
# gwas.regions <- list()
# gwas.regions[[1]] <- data.frame(gwas.hits[1, ])
#
# for(i in 2:nrow(GenOSept)){
#   print(i)
#   to.add <- TRUE
#     for(j in 1:nrow(gwas.loci)){
#         if(GenOSept$CHR[i] == gwas.loci$CHR[j] &
#            GenOSept$BP[i] %in% (gwas.loci$BP[j] - merging.distance):(gwas.loci$BP[j] + merging.distance)){
#             gwas.regions[[j]] <- rbind(gwas.regions[[j]], GenOSept[i, ])
#             to.add <- FALSE
#             break
#         }
#     }
#   if(to.add == TRUE & GenOSept$P.R.[i] < 1e-4) {
#     gwas.loci <- rbind(gwas.loci, GenOSept[i, ])
#     gwas.regions <- list.append(gwas.regions, GenOSept[i, ])
#   }
# }
#
# sum.eqtl.region <- function(x){
#   total <- sum(x$N.eQTL)
#   return(total)
# }
#
# gwas.loci$N.eQTL.region <- unlist(lapply(gwas.regions, sum.eqtl.region))
# to.test <- which(gwas.loci$N.eQTL.region > 0)
# gwas.loci <- gwas.loci[to.test, ]
# gwas.regions <- gwas.regions[to.test]
#
# write.table(gwas.loci, "../GWAS overlap/Peak_GWAS_hits.txt", sep="\t")

#################################################################################

# Overlap of Gene-SNP eqtl pairs
# Use data frame with

l.eqtl<-data.frame(fread("/Users/edaven/Documents/Data/Franke_eQTL_RNAseq/gene_level_eQTLs_z.txt"))
dim(l.eqtl) #5526971       4

l.eqtl$Gene.SNP<-paste(l.eqtl$ProbeName, l.eqtl$SNPName, sep="-")

tss.genes<-summ[,c("Ensembl_ID", "RsID", "eQTL_tvalue", "Minor_allele", "Major_allele")]

tss.genes$Gene.name<-do.call(rbind, strsplit(as.character(tss.genes$Ensembl_ID), "[.]"))[,1]
tss.genes$Gene.SNP<-paste(tss.genes$Gene.name, tss.genes$RsID, sep="-")

length(which(tss.genes$Gene.name%in%l.eqtl$ProbeName))
#4557 genes overlapping
length(which(tss.genes$Gene.SNP%in%l.eqtl$Gene.SNP))
#4250 gene SNP pairs overlapping = 85.4% overlap

eqtl.b<-merge(tss.genes, l.eqtl[, c("AlleleAssessed", "OverallZScore", "Gene.SNP")], by.x="Gene.SNP", by.y="Gene.SNP")
dim(eqtl.b) #4250 9

#Check which allele was tested in each study
eqtl.b$action<-"NA"
eqtl.b[which(eqtl.b$Minor_allele==eqtl.b$AlleleAssessed),"action"]<-"same"

table(eqtl.b$action)
#   NA same
# 2301 1949

#Allele we can't distinguish because of strand
eqtl.b[which(eqtl.b$Minor_allele=="A" & eqtl.b$Major_allele=="T" |
               eqtl.b$Minor_allele=="T" & eqtl.b$Major_allele=="A" |
               eqtl.b$Minor_allele=="C" & eqtl.b$Major_allele=="G" |
               eqtl.b$Minor_allele=="G" & eqtl.b$Major_allele=="C"), "action"] <-"strand"

table(eqtl.b$action)
#   NA   same strand
# 2270   1920     60

#Same but different strand

eqtl.b[which(eqtl.b$Minor_allele=="A" & eqtl.b$AlleleAssessed=="T" & eqtl.b$action!="strand" |
               eqtl.b$Minor_allele=="T" & eqtl.b$AlleleAssessed=="A" & eqtl.b$action!="strand" |
               eqtl.b$Minor_allele=="C" & eqtl.b$AlleleAssessed=="G" & eqtl.b$action!="strand" |
               eqtl.b$Minor_allele=="G" & eqtl.b$AlleleAssessed=="C" & eqtl.b$action!="strand"),"action"] <-"same"

table(eqtl.b$action)
#  NA   same strand
# 324   3866     60

eqtl.b[which(eqtl.b$action=="NA"),"action"] <-"flip"

table(eqtl.b$action)
# flip   same strand
#  324   3866     60

eqtl.b$t.value2<-eqtl.b$eQTL_tvalue
eqtl.b[which(eqtl.b$action=="flip"), "t.value2"]<-eqtl.b[which(eqtl.b$action=="flip"), "t.value2"]*(-1)

#Remove 60 SNPs which can't be distinguished because of strand
pdf("Overlap_lupus_eqtl_LudeRNAseq_strand_removed_tss.pdf", useDingbats=FALSE)
plot(eqtl.b[which(eqtl.b$action!="strand"), "t.value2"], eqtl.b[which(eqtl.b$action!="strand"), "OverallZScore"], xlab="Lupus z score", ylab="Zhernakova et al z score", pch=16, cex=0.6)
dev.off()

table(eqtl.b[which(eqtl.b$action!="strand"), "t.value2"]>0, eqtl.b[which(eqtl.b$action!="strand"), "OverallZScore"]>0)
#       FALSE TRUE
# FALSE  1956   18
# TRUE     18 2198

#4154/4190 = 99.1% same direction

cor.test(eqtl.b[which(eqtl.b$action!="strand"), "t.value2"], eqtl.b[which(eqtl.b$action!="strand"), "OverallZScore"])
#0.93
cor.test(eqtl.b[which(eqtl.b$action!="strand"), "t.value2"], eqtl.b[which(eqtl.b$action!="strand"), "OverallZScore"], method="spearman")

#################################################################################































genes.list <- list()
for(i in 1:22){
  genes.list[[i]] <- subset(genes, chrom == i)
}

#20253 3
#dataframe with gene names as row names and columns for chromosome number, start and stop position of gene
# for each chromosome

# dataframe with SNP names as row names and columns for chromosome number and SNP position
# for each chr
chr.bim.list <- list()
for(i in 1:22){
  chr.bim.list[[i]] <- data.frame(fread(paste("/well/jknight/Sepsis/Genotyping/Combined_Data/",
                                              i, ".merged.imputed.filtered.eqtl.bim", sep=""),
                                        header=FALSE))
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
}
