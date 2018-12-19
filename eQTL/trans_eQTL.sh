# Trans eQTL

# Prepare genotyping data following GTEX methods

~/plink --bfile AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed --extract AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed_noLD.prune.in --maf 0.05 --make-bed --out AD_736_multi_ethnic_chip_updated_trans_eQTL


#   --bfile AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed
#   --extract AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed_noLD.prune.in
#   --maf 0.05
#   --make-bed
#   --out AD_736_multi_ethnic_chip_updated_trans_eQTL
# 
# 64410 MB RAM detected; reserving 32205 MB for main workspace.
# 899740 variants loaded from .bim file.
# 470 people (180 males, 290 females) loaded from .fam.
# 470 phenotype values loaded from .fam.
# --extract: 281812 variants remaining.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 442 founders and 28 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 47489 het. haploid genotypes present (see
# AD_736_multi_ethnic_chip_updated_trans_eQTL.hh ); many commands treat these as
# missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.999198.
# 126180 variants removed due to minor allele threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 155632 variants and 470 people pass filters and QC.
# Among remaining phenotypes, 95 are cases and 375 are controls.
# --make-bed to AD_736_multi_ethnic_chip_updated_trans_eQTL.bed +
# AD_736_multi_ethnic_chip_updated_trans_eQTL.bim +
# AD_736_multi_ethnic_chip_updated_trans_eQTL.fam ... done.


sh update_build.sh AD_736_multi_ethnic_chip_updated_trans_eQTL Multi-EthnicGlobal_A1-b38.strand AD_736_multi_ethnic_chip_updated_trans_eQTL_b38

~/plink --bfile AD_736_multi_ethnic_chip_updated_trans_eQTL_b38 --recode A --out AD_736_multi_ethnic_chip_updated_trans_eQTL_b38

# Make Rdata file
load("eQTL.25PCs.RData")

library(data.table)
geno <- data.frame(fread("/well/jknight/AbuDhabiRNA/eQTL/Genotyping/AD_736_multi_ethnic_chip_updated_trans_eQTL_b38.raw"))
rownames(geno) <- geno[, 1]
geno[, 1:6] <- NULL
colnames(geno) <- gsub("X", "", colnames(geno))
colnames(geno) <- substr(colnames(geno), 1, nchar(colnames(geno))-2)
geno <- as.matrix(geno)
  
save(list=c("exprs", "geno", "info"), file = "/well/jknight/AbuDhabiRNA/eQTL/trans.eQTL.25PCs.RData")

#################################################################################
# R Script for mapping eQTL
#################################################################################

#!/usr/bin/env Rscript

# R script to map eQTL for the ith gene
# /apps/well/R/3.3.0/bin/Rscript --vanilla R.eQTL.r 1
args <- commandArgs(TRUE)
i <- as.integer(args[1])

# load R packages
suppressMessages(library(lme4))

# load necessary data
load("/well/jknight/AbuDhabiRNA/eQTL/trans.eQTL.25PCs.RData")
# exprs with 25 gex PCs and 4 geno PCs regressed out
# geno
# info

# output directory
outputdir <- "/well/jknight/AbuDhabiRNA/eQTL/Results/trans/"

# define the gene being tested
geneexprs <- exprs[i, ]
genename <- rownames(exprs)[i]

# pairs for this gene
irange <- 1:ncol(geno)

results <- do.call(rbind, lapply(irange, function(x){

        snp <- colnames(geno)[x]
        model <- lm(geneexprs ~ geno[, x])
        summary(model)$coefficients[2, ]

}))

colnames(results) <- c("eQTL_beta", "eQTL_SE", "eQTL_t", "eQTL_pval")

results <- data.frame(Gene=genename, SNP=colnames(geno), results)

saveRDS(results, paste0(outputdir, "result_final_gene", i, ".rds"))

#################################################################################
# Bash script for mapping eQTL
#################################################################################

#!/bin/bash

#$ -t 1:21580:40
#$ -cwd -V
#$ -N AD-trans
#$ -P jknight.prjc
#$ -q long.qc
#$ -o /well/jknight/AbuDhabiRNA/eQTL/Scripts/eQTL_trans_stdout
#$ -e /well/jknight/AbuDhabiRNA/eQTL/Scripts/eQTL_trans_stderr

echo "Starting the job $SGE_TASK_ID at" `date +"%T %m-%d-%Y"`
# Calculate the last task id for this step
this_step_last=$(( SGE_TASK_ID + SGE_TASK_STEPSIZE - 1 ))
if [ "${SGE_TASK_LAST}" -lt "${this_step_last}" ]
then
    this_step_last="${SGE_TASK_LAST}"
fi

# Loop over task ids in this step
while [ "${SGE_TASK_ID}" -le "${this_step_last}" ]
do
    echo `date +"%T %m-%d-%Y"`: SGE_TASK_ID=`printenv SGE_TASK_ID`
    /apps/well/R/3.3.0/bin/Rscript --vanilla  /well/jknight/AbuDhabiRNA/eQTL/Scripts/trans_eQTL.R $SGE_TASK_ID
    # Increment SGE_TASK_ID
    export SGE_TASK_ID=$(( SGE_TASK_ID + 1 ))
done
echo "************************"
echo "Finishing the job $this_step_last at" `date +"%T %m-%d-%Y"`
echo "************************"


#################################################################################
# R script for combining eQTL results and calculating FDR
#################################################################################

#!/usr/bin/env Rscript

cat.rds <- function(input_file_pattern, output_filename){
        temp_eqtl <- list.files(pattern = input_file_pattern)
        print(length(temp_eqtl))
        eqtl <- do.call(rbind, lapply(temp_eqtl, readRDS))
        print(dim(eqtl))
        saveRDS(eqtl, output_filename) 
}

cat.rds("/well/jknight/AbuDhabiRNA/eQTL/Results/trans/result_final_gene*", 
        "/well/jknight/AbuDhabiRNA/eQTL/Results/trans/eQTL_results_lm_trans.rds")

eqtl <- readRDS("/well/jknight/AbuDhabiRNA/eQTL/Results/trans/eQTL_results_lm_trans.rds")
eqtl$FDR <- p.adjust(eqtl$eQTL_pval, method="fdr")
saveRDS(eqtl.u, "/well/jknight/AbuDhabiRNA/eQTL/Results/trans/eQTL_results_lm_trans.rds")
