# R script to map interaction eQTL for eGenes using adjusted expression data
# /apps/well/R/3.3.0/bin/Rscript --vanilla R.eQTL.r "Disease"
mod <- "Disease"
# mod <- "SCORE.1.Emirati.0.Non"

# load R packages
# suppressMessages(library(lme4))

# load necessary data
load("/well/jknight/AbuDhabiRNA/eQTL/eQTL.25PCs.RData")
# exprs
# pairs
# geno
# info
sig.eQTL <- read.delim("/well/jknight/AbuDhabiRNA/eQTL/Results/cis/Peak_sig_eQTL_empP.txt", 
                       sep="\t", stringsAsFactors=F)
sig.eQTL <- readRDS("/well/jknight/AbuDhabiRNA/eQTL/Results/cis/Results_peak_final.rds")
sig.eQTL$Gene <- as.character(sig.eQTL$Gene)
sig.eQTL$SNP <- as.character(sig.eQTL$SNP)

# output directory
outputdir <- "/well/jknight/AbuDhabiRNA/eQTL/Results/cis/Interactions/"

# Test each gene
irange <- 1:nrow(sig.eQTL)
results <- do.call(rbind, lapply(irange, function(x){
  
  print(x)
  # define the gene being tested
  genename <- sig.eQTL$Gene[x]
  geneexprs <- exprs[match(genename, rownames(exprs)), ]

  # snp for this gene
  snp <- sig.eQTL$SNP[x]
  
  # model
  model <- lm(geneexprs ~ geno[, snp] + info[, mod] + geno[, snp]*info[, mod])
  
  if (length(which(geno[, snp]==2 & info[, mod] == "Control")) > 1 &
      #length(which(geno[, snp]==2 & info[, mod] == unique(info[, mod])[2])) > 1 &
      #length(which(geno[, snp]==2 & info[, mod] == unique(info[, mod])[3])) > 1 &
      length(which(geno[, snp]==2 & info[, mod] == "T2DM")) > 1){
  
  c(summary(model)$coefficients["geno[, snp]", ],
    summary(model)$coefficients["info[, mod]T2DM", ],
    summary(model)$coefficients["geno[, snp]:info[, mod]T2DM", ])
    
  } else {
    rep(NA, 12)
    
  }
  
}))

colnames(results)<-c("eQTL_beta", "eQTL_SE", "eQTL_t", "eQTL_pval",
                     "DE_beta", "DE_SE", "DE_t", "DE_pval",
                     "Interaction_beta", "Interaction_SE", "Interaction_t", "Interaction_pval")

results <- data.frame(Gene=sig.eQTL$Gene, SNP=sig.eQTL$SNP, results)

saveRDS(results, paste0(outputdir, "results_", mod, "_interactions.rds"))


################################################################################
################################################################################
################################################################################

mod <- "HbA1c"

results <- do.call(rbind, lapply(irange, function(x){
  
  print(x)
  # define the gene being tested
  genename <- sig.eQTL$Gene[x]
  geneexprs <- exprs[match(genename, rownames(exprs)), ]
  
  # snp for this gene
  snp <- sig.eQTL$SNP[x]
  
  # model
  model <- lm(geneexprs ~ geno[, snp] + info[, mod] + geno[, snp]*info[, mod])
  
  if ((length(which(geno[, snp]==0 & !(is.na(info[, mod])))) > 1) &
      (length(which(geno[, snp]==1 & !(is.na(info[, mod])))) > 1) &
      (length(which(geno[, snp]==2 & !(is.na(info[, mod])))) > 1)) {
    
    c(summary(model)$coefficients["geno[, snp]", ],
      summary(model)$coefficients["info[, mod]", ],
      summary(model)$coefficients["geno[, snp]:info[, mod]", ])
    
  } else {
    rep(NA, 12)
    
  }
  
}))

colnames(results)<-c("eQTL_beta", "eQTL_SE", "eQTL_t", "eQTL_pval",
                     "DE_beta", "DE_SE", "DE_t", "DE_pval",
                     "Interaction_beta", "Interaction_SE", "Interaction_t", "Interaction_pval")

results <- data.frame(Gene=sig.eQTL$Gene, SNP=sig.eQTL$SNP, results)

saveRDS(results, paste0(outputdir, "results_", mod, "_interactions.rds"))



################################################################################
################################################################################
################################################################################


mod <- "Glucose"

results <- do.call(rbind, lapply(irange, function(x){
  
  print(x)
  # define the gene being tested
  genename <- sig.eQTL$Gene[x]
  geneexprs <- exprs[match(genename, rownames(exprs)), ]
  
  # snp for this gene
  snp <- sig.eQTL$SNP[x]
  
  # model
  model <- lm(geneexprs ~ geno[, snp] + info[, mod] + geno[, snp]*info[, mod])
  
  if ((length(which(geno[, snp]==0 & !(is.na(info[, mod])))) > 1) &
    (length(which(geno[, snp]==1 & !(is.na(info[, mod])))) > 1) &
  (length(which(geno[, snp]==2 & !(is.na(info[, mod])))) > 1)) {
    
    c(summary(model)$coefficients["geno[, snp]", ],
      summary(model)$coefficients["info[, mod]", ],
      summary(model)$coefficients["geno[, snp]:info[, mod]", ])
    
  } else {
    rep(NA, 12)
    
  }
  
}))

colnames(results)<-c("eQTL_beta", "eQTL_SE", "eQTL_t", "eQTL_pval",
                     "DE_beta", "DE_SE", "DE_t", "DE_pval",
                     "Interaction_beta", "Interaction_SE", "Interaction_t", "Interaction_pval")

results <- data.frame(Gene=sig.eQTL$Gene, SNP=sig.eQTL$SNP, results)

saveRDS(results, paste0(outputdir, "results_", mod, "_interactions.rds"))






