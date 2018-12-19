# Identify significant trans eQTL

n.trans <- 21580 * 155622
n.cis <- 3473136

p.val <- 0.05/(n.trans - n.cis)

input_file_pattern <- "result_final_gene*"
temp_eqtl <- list.files(pattern = input_file_pattern)
print(length(temp_eqtl))

cis.eqtl <- readRDS("/well/jknight/AbuDhabiRNA/eQTL/Results/cis/eQTL_results_final_lm_cis.rds")
cis.eqtl$Gene <- as.character(cis.eqtl$Gene)
cis.eqtl$SNP <- as.character(cis.eqtl$SNP)
cis.eqtl$GeneSNP <- paste0(cis.eqtl$Gene, cis.eqtl$SNP)

trans.eqtl <- lapply(temp_eqtl, function(x){
  print(x)
  gene <- readRDS(x)
  print(dim(gene))
  gene$GeneSNP <- paste0(gene$Gene, gene$SNP)
  gene <- gene[!(gene$GeneSNP %in% cis.eqtl$GeneSNP), ]
  print(dim(gene))
  sig <- gene[gene$eQTL_pval <= p.val, ]
  print(dim(sig))
  return(sig)
})

sig.eqtl <- do.call(rbind, trans.eqtl)
print(dim(sig.eqtl))

saveRDS(eqtl, "/well/jknight/AbuDhabiRNA/eQTL/Results/trans/eQTL_results_lm_sig_trans.rds") 
