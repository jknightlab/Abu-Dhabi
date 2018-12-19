# Population-specific SNPs MAF

library(data.table)

AD <- fread("/well/jknight/AbuDhabiRNA/eQTL/Genotyping/AD2463.frq")
bim <- fread("/well/jknight/AbuDhabiRNA/eQTL/Genotyping/AD2463.bim")
AD$pos.id <- paste0(AD$CHR, "-", bim$V4[match(AD$SNP, bim$V2)])

# Columns are pos-id, rsid, ref, alt, ref_count, alt_count.

gnomad.files <- list.files(path="/well/jknight/AbuDhabiRNA/eQTL/Genotyping/gnomAD",
                     pattern="gnomad_AFR*")

gnomad.afr <- lapply(gnomad.files, fread)

gnomad.files <- list.files(path="/well/jknight/AbuDhabiRNA/eQTL/Genotyping/gnomAD",
                           pattern="gnomad_EAS*")

gnomad.eas <- lapply(gnomad.files, fread)

gnomad.files <- list.files(path="/well/jknight/AbuDhabiRNA/eQTL/Genotyping/gnomAD",
                           pattern="gnomad_NFE*")

gnomad.nfe <- lapply(gnomad.files, fread)

#################################
# Match
############

gnomad.afr <- lapply(gnomad.afr, function(x){
  out <- x[x$V1 %in% AD$pos.id]
  return(out)
})
gnomad.afr <- rbindlist(gnomad.afr)

gnomad.eas <- lapply(gnomad.eas, function(x){
  out <- x[x$V1 %in% AD$pos.id]
  return(out)
})
gnomad.eas <- rbindlist(gnomad.eas)

gnomad.nfe <- lapply(gnomad.nfe, function(x){
  out <- x[x$V1 %in% AD$pos.id]
  return(out)
})
gnomad.nfe <- rbindlist(gnomad.nfe)

AD <- AD[AD$pos.id %in% gnomad.afr$V1, ]
gnomad.afr <- gnomad.afr[match(AD$pos.id, gnomad.afr$V1), ]
gnomad.eas <- gnomad.eas[match(AD$pos.id, gnomad.eas$V1), ]
gnomad.nfe <- gnomad.nfe[match(AD$pos.id, gnomad.nfe$V1), ]

gnomad.afr$AFR <- gnomad.afr$V6/(gnomad.afr$V5 + gnomad.afr$V6)
gnomad.eas$EAS <- gnomad.eas$V6/(gnomad.eas$V5 + gnomad.eas$V6)
gnomad.nfe$NFE <- gnomad.nfe$V6/(gnomad.nfe$V5 + gnomad.nfe$V6)

# check alleles
gnomad.afr$action <- "NA"
gnomad.afr[which(gnomad.afr$V4 == AD$A1), "action"] <- "same"

table(gnomad.afr$action)
# NA  same 
# 18341 17090 

# Alleles can't distinguish because of strand
gnomad.afr[which(gnomad.afr$V4=="A" & gnomad.afr$V3=="T" |
                     gnomad.afr$V4=="T" & gnomad.afr$V3=="A" |
                     gnomad.afr$V4=="C" & gnomad.afr$V3=="G" |
                     gnomad.afr$V4=="G" & gnomad.afr$V3=="C"), 
             "action"] <-"strand"

table(gnomad.afr$action)
# NA   same strand 
# 18198  16962    271 

# Same but different strand
gnomad.afr[which(gnomad.afr$V4=="A" & AD$A1=="T" & 
                     gnomad.afr$action != "strand" |
                     gnomad.afr$V4=="T" & AD$A1=="A" & 
                     gnomad.afr$action != "strand" |
                     gnomad.afr$V4=="C" & AD$A1=="G" & 
                     gnomad.afr$action != "strand" |
                     gnomad.afr$V4=="G" & AD$A1=="C" & 
                     gnomad.afr$action != "strand"),
             "action"] <-"same"

table(gnomad.afr$action)
# NA   same strand 
# 1425  33735    271  

gnomad.afr[which(gnomad.afr$action == "NA"), "action"] <- "flip"

table(gnomad.afr$action)
# flip   same strand 
# 1425  33735    271 

gnomad.afr$maf <- gnomad.afr$AFR
gnomad.afr[which(gnomad.afr$action == "flip"), "maf"] <- (1 - 
  gnomad.afr[which(gnomad.afr$action == "flip"), "maf"])

# Remove SNPs which can't be distinguished because of strand
gnomad.afr$maf[which(gnomad.afr$action == "strand")] <- NA
gnomad.afr <- gnomad.afr[complete.cases(gnomad.afr), ]




gnomad.eas$action <- "NA"
gnomad.eas[which(gnomad.eas$V4 == AD$A1), "action"] <- "same"

table(gnomad.eas$action)
# NA  same 
# 18341 17090 

# Alleles can't distinguish because of strand
gnomad.eas[which(gnomad.eas$V4=="A" & gnomad.eas$V3=="T" |
                   gnomad.eas$V4=="T" & gnomad.eas$V3=="A" |
                   gnomad.eas$V4=="C" & gnomad.eas$V3=="G" |
                   gnomad.eas$V4=="G" & gnomad.eas$V3=="C"), 
           "action"] <-"strand"

table(gnomad.eas$action)
# NA   same strand 
# 18198  16962    271 

# Same but different strand
gnomad.eas[which(gnomad.eas$V4=="A" & AD$A1=="T" & 
                   gnomad.eas$action != "strand" |
                   gnomad.eas$V4=="T" & AD$A1=="A" & 
                   gnomad.eas$action != "strand" |
                   gnomad.eas$V4=="C" & AD$A1=="G" & 
                   gnomad.eas$action != "strand" |
                   gnomad.eas$V4=="G" & AD$A1=="C" & 
                   gnomad.eas$action != "strand"),
           "action"] <-"same"

table(gnomad.eas$action)
# NA   same strand 
# 1425  33735    271  

gnomad.eas[which(gnomad.eas$action == "NA"), "action"] <- "flip"

table(gnomad.eas$action)
# flip   same strand 
# 1425  33735    271 

gnomad.eas$maf <- gnomad.eas$EAS
gnomad.eas[which(gnomad.eas$action == "flip"), "maf"] <- (1 - 
                                                            gnomad.eas[which(gnomad.eas$action == "flip"), "maf"])

# Remove SNPs which can't be distinguished because of strand
gnomad.eas$maf[which(gnomad.eas$action == "strand")] <- NA
gnomad.eas <- gnomad.eas[complete.cases(gnomad.eas), ]





gnomad.nfe$action <- "NA"
gnomad.nfe[which(gnomad.nfe$V4 == AD$A1), "action"] <- "same"

table(gnomad.nfe$action)
# NA  same 
# 18341 17090 

# Alleles can't distinguish because of strand
gnomad.nfe[which(gnomad.nfe$V4=="A" & gnomad.nfe$V3=="T" |
                   gnomad.nfe$V4=="T" & gnomad.nfe$V3=="A" |
                   gnomad.nfe$V4=="C" & gnomad.nfe$V3=="G" |
                   gnomad.nfe$V4=="G" & gnomad.nfe$V3=="C"), 
           "action"] <-"strand"

table(gnomad.nfe$action)
# NA   same strand 
# 18198  16962    271 

# Same but different strand
gnomad.nfe[which(gnomad.nfe$V4=="A" & AD$A1=="T" & 
                   gnomad.nfe$action != "strand" |
                   gnomad.nfe$V4=="T" & AD$A1=="A" & 
                   gnomad.nfe$action != "strand" |
                   gnomad.nfe$V4=="C" & AD$A1=="G" & 
                   gnomad.nfe$action != "strand" |
                   gnomad.nfe$V4=="G" & AD$A1=="C" & 
                   gnomad.nfe$action != "strand"),
           "action"] <-"same"

table(gnomad.nfe$action)
# NA   same strand 
# 1425  33735    271  

gnomad.nfe[which(gnomad.nfe$action == "NA"), "action"] <- "flip"

table(gnomad.nfe$action)
# flip   same strand 
# 1425  33735    271 

gnomad.nfe$maf <- gnomad.nfe$NFE
gnomad.nfe[which(gnomad.nfe$action == "flip"), "maf"] <- (1 - 
                                                            gnomad.nfe[which(gnomad.nfe$action == "flip"), "maf"])

# Remove SNPs which can't be distinguished because of strand
gnomad.nfe$maf[which(gnomad.nfe$action == "strand")] <- NA
gnomad.nfe <- gnomad.nfe[complete.cases(gnomad.nfe), ]


AD$AFR <- gnomad.afr$maf[match(AD$pos.id, gnomad.afr$V1)]
AD$EAS <- gnomad.eas$maf[match(AD$pos.id, gnomad.eas$V1)]
AD$NFE <- gnomad.nfe$maf[match(AD$pos.id, gnomad.nfe$V1)]

#################################################################################

maf <- read.delim("../eQTL/MAF_Comparison.txt")
load("../eQTL/ShinyApp/Peak_eQTL.RData")
maf <- maf[maf$SNP %in% peak.eQTL$SNP[peak.eQTL$qvalue < 0.05], ]

ggplot(maf, aes(MAF, AFR)) + 
  theme_bw() +
   geom_point() +
   geom_abline(intercept=0, slope=1)

ggplot(maf, aes(MAF, EAS)) + 
  theme_bw() +
  geom_point() +
  geom_abline(intercept=0, slope=1)

ggplot(maf, aes(MAF, NFE)) + 
  theme_bw() +
  geom_point() +
  geom_abline(intercept=0, slope=1)

ggplot(maf, aes(EAS, AFR)) + 
  theme_bw() +
  geom_point() +
  geom_abline(intercept=0, slope=1)

ggplot(maf, aes(EAS, NFE)) + 
  theme_bw() +
  geom_point() +
  geom_abline(intercept=0, slope=1)

all.eqtl$MAF <- maf$MAF[match(all.eqtl$SNP, maf$SNP)]

ggplot(all.eqtl, aes(MAF, abs(eQTL_beta))) +
  geom_smooth() +
   theme_bw()
