# QC of Abu Dhabi genotyping data

Take as input files from Kate: AD_736_multi_ethnic_chip_update_strand

Update family info and IDs so that they match the RNA-seq data:

```{bash}
~/plink --bfile AD_736_multi_ethnic_chip_update_strand --update-parents update_parents.txt --make-bed --out AD_736_multi_ethnic_chip_update_strand_update_parents

~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents --update-ids update_FID_info.txt --make-bed --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID
```

# Sex check

There are some incorrect sex codes - update-sex.txt prepared based on the most
recent clinical data.

```{bash}
~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID --update-sex update-sex.txt --make-bed --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex
```

From the Plink manual the following steps are recommended for checking the sample's
recorded sex against the sex predicted from the genotyping data:

```{bash}
~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex --split-x hg19 --make-bed --out AD_736_multi_ethnic_chip_updated_for_sex_check

~/plink --bfile AD_736_multi_ethnic_chip_updated_for_sex_check --indep-pairphase 20000 2000 0.5 --out AD_736_multi_ethnic_chip_updated_for_sex_check

~/plink --bfile AD_736_multi_ethnic_chip_updated_for_sex_check --extract AD_736_multi_ethnic_chip_updated_for_sex_check.prune.in --check-sex --out AD_736_multi_ethnic_chip_updated_sex_check

grep PROBLEM AD_736_multi_ethnic_chip_updated_sex_check.sexcheck > sexprobs

awk '{print $1}' sexprobs > samples_to_remove_sex.txt
```

There are 31 problems:
       FG83         F83            1            2      PROBLEM       0.1072
       FG84         F84            1            2      PROBLEM     -0.03324
       FG13        FG13            2            0      PROBLEM       0.2477
       FG55        FG55            2            0      PROBLEM        0.263
       FG71        FG71            2            0      PROBLEM        0.263
      FG151       FG151            2            0      PROBLEM       0.4014
       FG39         M39            2            1      PROBLEM       0.9446
       FG84         M84            2            1      PROBLEM       0.9492
        S70         S70            1            0      PROBLEM       0.5563
       S116        S116            2            0      PROBLEM       0.2776
       S124        S124            2            0      PROBLEM        0.348
       S157        S157            2            0      PROBLEM       0.2038
       S168        S168            2            0      PROBLEM       0.2015
       S185        S185            2            0      PROBLEM       0.2917
       S202        S202            2            0      PROBLEM       0.2155
       S221        S221            2            0      PROBLEM       0.3436
       S222        S222            2            0      PROBLEM       0.3383
       S227        S227            2            0      PROBLEM       0.5283
       S278        S278            2            0      PROBLEM       0.2232
       S308        S308            2            0      PROBLEM       0.3258
       S313        S313            2            0      PROBLEM       0.2538
       S315        S315            2            0      PROBLEM       0.2125
       S363        S363            2            0      PROBLEM       0.3251
       S395        S395            2            1      PROBLEM       0.9474
       S414        S414            2            0      PROBLEM        0.335
       S416        S416            2            0      PROBLEM       0.2094
       S419        S419            2            0      PROBLEM       0.2029
       S440        S440            2            0      PROBLEM       0.2142
       S449        S449            2            0      PROBLEM       0.2418
       S458        S458            2            0      PROBLEM       0.3255
       S480        S480            1            2      PROBLEM       0.1029

Plot the results in R and redefine the cut-offs based on the data:

```{r}
sex.check <- read.delim("AD_736_multi_ethnic_chip_update_strand_sex_check.sexcheck", sep="")
pdf("Sex_check_F_estimates.pdf")
hist(sex.check$F, breaks=20)
dev.off()
```

```{bash}
~/plink --bfile AD_736_multi_ethnic_chip_updated_for_sex_check --extract AD_736_multi_ethnic_chip_updated_for_sex_check.prune.in --check-sex 0.6 0.8 --out AD_736_multi_ethnic_chip_updated_sex_check

grep PROBLEM AD_736_multi_ethnic_chip_updated_sex_check.sexcheck > sexprobs

awk '{print $1, $2}' sexprobs > fail_sex.txt
```

7 problems:

FID				ID		PEDSEX	SNPSEX	STATUS	F_ESTIMATE					NOTES
FG83         F83            1            2      PROBLEM       0.1072			
FG84         F84            1            2      PROBLEM     -0.03324			
FG39         M39            2            1      PROBLEM       0.9446			
FG84         M84            2            1      PROBLEM       0.9492			
S70         S70            1            2      PROBLEM       0.5563				In sample information file male - exclude
S395        S395            2            1      PROBLEM       0.9474			In sample information file female - EXCLUDE
S480        S480            1            2      PROBLEM       0.1029			In sample information file male - EXCLUDE


FG83, FG84, FG39, FG84 are parents labelled the wrong way round - update labels, sex
 
Update family info:

```{bash}
~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex --update-ids update_FID_info_postsexcheck.txt --make-bed --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked

~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked --update-sex update-sex_postsexcheck.txt --make-bed --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked
```

S395 and S480 are also misgendered in the GEx - worth getting Hinda to check?
For now excluded with S70.

```{bash}
~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked --split-x hg19 --make-bed --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID__for_sex_check

~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID__for_sex_check --extract AD_736_multi_ethnic_chip_updated_for_sex_check.prune.in --check-sex 0.6 0.7 --out AD_736_multi_ethnic_chip_updated_sex_check

grep PROBLEM AD_736_multi_ethnic_chip_updated_sex_check.sexcheck > sexprobs
awk '{print $1, $2}' sexprobs > fail_sex.txt
```

# Missingness

```{bash}
~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked --missing --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked

~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked --het --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked
```

```{r}
imiss <- read.table("AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked.imiss", h=T)
imiss$logF_MISS <- log10(imiss[, 6])

het <- read.table("AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked.het", h=T)
het$meanHet <- (het$N.NM. - het$O.HOM.) / het$N.NM.

colors  <- densCols(imiss$logF_MISS, het$meanHet)

y.min <- min(c(het$meanHet, (mean(het$meanHet) - (3 * sd(het$meanHet)))))
y.max <- max(c(het$meanHet, (mean(het$meanHet) + (3 * sd(het$meanHet)))))

pdf("AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_missing_het.pdf")
plot(imiss$logF_MISS, het$meanHet, col=colors, pch=20, xlab="Proportion of missing genotypes", ylab="Heterozygosity rate", ylim=c(y.min, y.max), axes=F)
axis(1, at=c(-3,-2,-1,0), labels=c(0.001,0.01,0.1,1))
axis(2, at=c(0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), tick=T)
abline(h=mean(het$meanHet) - (3 * sd(het$meanHet)), col="RED", lty=2)
abline(h=mean(het$meanHet) + (3 * sd(het$meanHet)), col="RED", lty=2)
abline(v=log10(0.03), col="RED", lty=2)
dev.off()

missing <- imiss[which(imiss$logF_MISS > log10(0.03)), 1:2]
print(missing)
# "FG64" "FG69" "S70"  "S97"  "S191"

write.table(missing, "fail_missing", sep="\t", row.names=F, quote=F, col.names=F)

bad.het <- het[which(het$meanHet < mean(het$meanHet) - (3 * sd(het$meanHet)) |
        het$meanHet > mean(het$meanHet) + (3 * sd(het$meanHet))), 1:2]
print(bad.het)
# S191
write.table(bad.het, "fail_het", sep="\t", row.names=F, quote=F, col.names=F)
```

# Relatedness

This is complicated due to the population being studied. We have therefore 
decided to exclude parents and other similarly close relatives i.e. pi_hat >=0.5.

```{bash}
~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked --extract AD_736_multi_ethnic_chip_updated_for_sex_check.prune.in --genome --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked
```

```{r}
ibd <- read.delim("AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked.genome", sep="")
pi_hat <- subset(ibd, PI_HAT > 0.5) #84

# Lots of problems - exclude one from each pair with >0.5
# This also includes parents

fam <- read.delim("AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked.fam", 
sep="", header=F, stringsAsFactors=F)

parents <- which(fam$V1 != fam$V2)
to.rm <- fam[parents, 1:2] #86

pi_hat <- pi_hat[!(pi_hat$IID1 %in% to.rm$V2 | pi_hat$IID2 %in% to.rm$V2), ] # 14
outliers <- unique(c("FG64", "FG69", "S70",  "S97",  "S191", "S70", "S395", "S480"))
pi_hat <- pi_hat[!(pi_hat$IID1 %in% outliers | pi_hat$IID2 %in% outliers), ]

write.table(pi_hat, "pi_hat.txt", sep="\t", quote=F)
colnames(to.rm) <- c("FID1", "IID1")
write.table(to.rm, "Just_Parents_to_rm.txt", sep="\t", row.names=FALSE, quote=F, col.names=F)
to.rm <- rbind(to.rm, pi_hat[, 1:2])
write.table(to.rm, "Parents_to_rm.txt", sep="\t", row.names=FALSE, quote=F, col.names=F)
```

14 pairs:

        FID1  IID1 FID2 IID2 RT EZ     Z0     Z1     Z2 PI_HAT PHE      DST PPC
59406   FG52  FG52 S352 S352 UN NA 0.1774 0.4809 0.3417 0.5821   0 0.941484   1
116660 FG147 FG147 S416 S416 UN NA 0.2176 0.4299 0.3525 0.5674  -1 0.940292   1
149830    S8    S8 S146 S146 UN NA 0.1278 0.4681 0.4041 0.6381   0 0.948664   1
168753   S50   S50 S248 S248 UN NA 0.2307 0.4679 0.3014 0.5353   0 0.935784   1
196176  S116  S116 S117 S117 UN NA 0.0070 0.9820 0.0110 0.5020   1 0.925039   1
197710  S120  S120 S121 S121 UN NA 0.0041 0.8918 0.1040 0.5499   0 0.932199   1
209921  S153  S153 S320 S320 UN NA 0.2106 0.5238 0.2656 0.5275  -1 0.934084   1
221340  S188  S188 S189 S189 UN NA 0.2204 0.4409 0.3388 0.5592  -1 0.939122   1
227904  S208  S208 S411 S411 UN NA 0.2021 0.4963 0.3015 0.5497  -1 0.937220   1
255081  S326  S326 S327 S327 UN NA 0.0123 0.9300 0.0577 0.5227  -1 0.928296   1
256348  S333  S333 S397 S397 UN NA 0.0127 0.9625 0.0248 0.5061  -1 0.925796   1
256453  S334  S334 S335 S335 UN NA 0.1870 0.4562 0.3568 0.5849  -1 0.942146   1
265228  S399  S399 S400 S400 UN NA 0.0083 0.9281 0.0636 0.5276  -1 0.928938   1
267003  S418  S418 S427 S427 UN NA 0.2374 0.5160 0.2466 0.5046   0 0.931314   1

These generally look like they could be siblings? They are frequently recruited
at the same time/have consecutive sample numbers.

At this stage, cross reference with the gene expression data to see if any have
already been excluded there:

-> S326 (IID1)

I will run the PCA, see where these individuals are - probably looking weird in pairs

# MDS

```{bash}
~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked --read-genome AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked.genome --cluster --mds-plot 4 --extract AD_736_multi_ethnic_chip_updated_for_sex_check.prune.in --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked
```

```{r}
mds <- read.delim("AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked.mds", sep="")
mds$QC_removal <- "None"

sex <- read.delim("fail_sex.txt", header=F)
mds$QC_removal[match(sex$V1, mds$IID)] <- "Sex"

missing <- read.delim("fail_missing", header=F)
mds$QC_removal[match(missing$V1, mds$IID)] <- "Missingness"

related <- read.delim("Parents_to_rm.txt")
mds$QC_removal[match(related$IID, mds$IID)] <- "Relatedness"

mds$QC_removal <- as.factor(mds$QC_removal)

pdf("mds_plots_all.pdf", onefile=T, width=11, height=9)

palette(c("red", "blue", "black", "green"))
plot(mds$C1, mds$C2, col=mds$QC_removal, xlab="First dimension", 
     ylab="Second dimension", pch=16)
legend("topright", legend=unique(mds$QC_removal), 
       pch=16, col=unique(mds$QC_removal), cex=0.7, 
       title="QC removal reason")

plot(mds$C1, mds$C3, col=mds$QC_removal, xlab="First dimension", 
     ylab="Third dimension", pch=16)
legend("topright", legend=unique(mds$QC_removal), 
       pch=16, col=unique(mds$QC_removal), cex=0.7, title="QC removal reason")

plot(mds$C2, mds$C3, col=mds$QC_removal, xlab="Second dimension", 
     ylab="Third dimension", pch=16)
legend("topright", legend=unique(mds$QC_removal), 
       pch=16, col=unique(mds$QC_removal), cex=0.7, title="QC removal reason")

plot(mds$C3, mds$C4, col=mds$QC_removal, xlab="Third dimension", 
     ylab="Fourth dimension", pch=16)
legend("topright", legend=unique(mds$QC_removal), 
       pch=16, col=unique(mds$QC_removal), cex=0.7, title="QC removal reason")

mds <- subset(mds, QC_removal == "None")
mds <- droplevels(mds)

plot(mds$C1, mds$C2, col=mds$QC_removal, xlab="First dimension", 
     ylab="Second dimension", pch=16)

plot(mds$C1, mds$C3, col=mds$QC_removal, xlab="First dimension", 
     ylab="Third dimension", pch=16)

plot(mds$C2, mds$C3, col=mds$QC_removal, xlab="Second dimension", 
     ylab="Third dimension", pch=16)

plot(mds$C3, mds$C4, col=mds$QC_removal, xlab="Third dimension", 
     ylab="Fourth dimension", pch=16)
dev.off()
```

```{bash}
# Without excluding pi_hat
cat fail_het fail_missing fail_sex.txt | sort -k1 | uniq > fail-qc-inds.txt
# 8
cat fail-qc-inds.txt Just_Parents_to_rm.txt | sort -k1 | uniq > fail-qc-inds2.txt
# 94 inds

# MDS
~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked --remove fail-qc-inds2.txt --extract AD_736_multi_ethnic_chip_updated_for_sex_check.prune.in --genome --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_first_ind_removed

~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked --read-genome AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_first_ind_removed.genome --remove fail-qc-inds2.txt --cluster --mds-plot 10 --extract AD_736_multi_ethnic_chip_updated_for_sex_check.prune.in --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_no_fails


#PCA
~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked --remove fail-qc-inds2.txt --exclude /well/jknight/Sepsis/Genotyping/high_LD_regions_hg19.txt --range --make-bed --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_noLongRangeLD

~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_noLongRangeLD --extract AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed_noLD.prune.in --pca --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_noLongRangeLD

```

PCA looks like it identifies the same individuals.

```{r}
geno.pc <- read.table("AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_noLongRangeLD.eigenvec", sep="", header=F, row.names=2)
geno.pc$V1 <- NULL
plot(geno.pc$V3, geno.pc$V4, pch=16)
plot(geno.pc$V5, geno.pc$V6, pch=16)
plot(geno.pc$V7, geno.pc$V8, pch=16)
rownames(geno.pc)[identify(geno.pc$V7, geno.pc$V8)]



mds <- read.delim("AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_no_fails.mds", sep="")
mds$QC_removal <- "None"
pi_hat <- read.delim("pi_hat.txt", header=T)
mds$QC_removal[match(pi_hat$IID1, mds$IID)] <- "IBD"
mds$QC_removal[match(pi_hat$IID2, mds$IID)] <- "IBD"

mds$QC_removal <- as.factor(mds$QC_removal)

pdf("mds_plots_no_fails.pdf", onefile=T, width=11, height=9)

plot(mds$C1, mds$C2, col=mds$QC_removal, xlab="First dimension", 
     ylab="Second dimension", pch=16)
plot(mds$C1, mds$C3, col=mds$QC_removal, xlab="First dimension", 
     ylab="Third dimension", pch=16)
plot(mds$C2, mds$C3, col=mds$QC_removal, xlab="Second dimension", 
     ylab="Third dimension", pch=16)
plot(mds$C3, mds$C4, col=mds$QC_removal, xlab="Third dimension", 
     ylab="Fourth dimension", pch=16)
 
mds[mds$C4 < -0.025, ]

      FID   IID SOL        C1         C2         C3         C4         C5         C6          C7         C8           C9
101 FG106 FG106   0 0.0383554 0.00529019 0.00673513 -0.0302652 -0.0809358 0.00469959 0.000629565 0.00143733 -0.001888040
450  S306  S306   0 0.0375259 0.00504339 0.00601879 -0.0305285 -0.0816016 0.00561861 0.001657870 0.00160944 -0.001576990
581  S439  S439   0 0.0375076 0.00505722 0.00535018 -0.0294851 -0.0786758 0.00677340 0.001290870 0.00166526 -0.000946804
            C10 QC_removal
101 -0.00347146       None
450 -0.00299845       None
581 -0.00161757       None

plot(mds$C5, mds$C6, col=mds$QC_removal, xlab="Fifth dimension", 
     ylab="Sixth dimension", pch=16)

mds[identify(mds$C5, mds$C6, labels=mds$IID), ]

      FID   IID SOL        C1         C2         C3          C4          C5          C6           C7         C8           C9
101 FG106 FG106   0 0.0383554 0.00529019 0.00673513 -0.03026520 -0.08093580  0.00469959  0.000629565 0.00143733 -0.001888040
202   S55   S55   0 0.0659774 0.01494970 0.00603186 -0.00253807  0.01392020 -0.00237898  0.045843800 0.03562450 -0.047632200
320  S175  S175   0 0.0679081 0.01549780 0.00570497 -0.00304744  0.01404330 -0.00213555  0.045818100 0.03444620 -0.047681000
359  S214  S214   0 0.0409481 0.01366690 0.00213783  0.00269103  0.00515379  0.00939991 -0.057797700 0.04896200  0.029288000
360  S215  S215   0 0.0421773 0.01354790 0.00211610  0.00325240  0.00572021  0.00918486 -0.057678300 0.04873950  0.030313500
450  S306  S306   0 0.0375259 0.00504339 0.00601879 -0.03052850 -0.08160160  0.00561861  0.001657870 0.00160944 -0.001576990
581  S439  S439   0 0.0375076 0.00505722 0.00535018 -0.02948510 -0.07867580  0.00677340  0.001290870 0.00166526 -0.000946804
            C10 QC_removal
101 -0.00347146       None
202 -0.04166380       None
320 -0.04140890       None
359 -0.03787980       None
360 -0.03802550       None
450 -0.00299845       None
581 -0.00161757       None

plot(mds$C7, mds$C8, col=mds$QC_removal, xlab="Seventh dimension", 
     ylab="Eighth dimension", pch=16)

mds[identify(mds$C7, mds$C8, labels=mds$IID), ]

     FID  IID SOL        C1          C2           C3          C4          C5          C6         C7          C8           C9
11  FG11 FG11   0 0.0103359 -0.02306940 -0.015997800  0.00103905 0.000727652  0.00216467  0.0217981  0.00189043  0.030696500
177  S29  S29   0 0.0554066  0.01362070  0.003770680 -0.00290161 0.007746530 -0.00149691 -0.0175880 -0.06902650  0.000338973
178  S30  S30   0 0.0528803  0.01306710  0.003539120 -0.00313662 0.008272440 -0.00113805 -0.0166575 -0.06946110  0.000361373
202  S55  S55   0 0.0659774  0.01494970  0.006031860 -0.00253807 0.013920200 -0.00237898  0.0458438  0.03562450 -0.047632200
320 S175 S175   0 0.0679081  0.01549780  0.005704970 -0.00304744 0.014043300 -0.00213555  0.0458181  0.03444620 -0.047681000
359 S214 S214   0 0.0409481  0.01366690  0.002137830  0.00269103 0.005153790  0.00939991 -0.0577977  0.04896200  0.029288000
360 S215 S215   0 0.0421773  0.01354790  0.002116100  0.00325240 0.005720210  0.00918486 -0.0576783  0.04873950  0.030313500
494 S350 S350   0 0.0133378  0.00964422  0.000335183  0.00766197 0.003198810  0.00337534 -0.0180694  0.00648067 -0.005923170
495 S351 S351   0 0.0403617  0.01253930  0.002085440  0.00366115 0.005883160  0.00287251 -0.0178353  0.00455386 -0.005076170
            C10 QC_removal
11  -0.00706772       None
177 -0.04943010       None
178 -0.04940130       None
202 -0.04166380       None
320 -0.04140890       None
359 -0.03787980       None
360 -0.03802550       None
494  0.02205140       None
495  0.02628330       None

plot(mds$C9, mds$C10, col=mds$QC_removal, xlab="Ninth dimension", 
     ylab="Tenth dimension", pch=16)
     
mds[identify(mds$C9, mds$C10, labels=mds$IID), ]

     FID  IID SOL        C1          C2           C3           C4          C5          C6         C7           C8          C9
11  FG11 FG11   0 0.0103025 -0.02364380 -0.015648500  0.000970638 0.000616365  0.00335050 -0.0190665  0.000280708 -0.03310950
174  S29  S29   0 0.0553443  0.01377090  0.003314680 -0.002572340 0.007797610  0.00157821  0.0192995  0.068826400  0.00727036
198  S55  S55   0 0.0659108  0.01518190  0.005523510 -0.002160740 0.014093600  0.01515300 -0.0467782 -0.038689400  0.04268600
313 S175 S175   0 0.0678426  0.01572220  0.005224900 -0.002635480 0.014211500  0.01532370 -0.0467376 -0.037534000  0.04261440
350 S214 S214   0 0.0408649  0.01374170  0.001669090  0.002982890 0.004936990 -0.03183100  0.0535058 -0.046094600 -0.02613780
482 S350 S350   0 0.0132485  0.00954125  0.000206884  0.007802540 0.002847410 -0.00929190  0.0158929 -0.006526010  0.00917633
483 S351 S351   0 0.0402828  0.01260200  0.001670620  0.003951040 0.005655830 -0.00896619  0.0157545 -0.004553960  0.00836644
            C10 QC_removal
11   0.00821355       None
174  0.04919020       None
198  0.04177560       None
313  0.04133100       None
350  0.03844960       None
482 -0.02374440       None
483 -0.02802640       None

dev.off()
```

# With pi_hat excluded
# 108 inds
cat fail-qc-inds.txt Parents_to_rm.txt | sort -k1 | uniq > fail-qc-inds2.txt
fail-qc-inds2.txt

# Repeat having exluded Pi_hat

mds <- read.delim("AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_no_fails_no_IBD.mds", sep="")
mds$QC_removal <- "None"
pi_hat <- read.delim("pi_hat.txt", header=T)
mds$QC_removal[match(pi_hat$IID1, mds$IID)] <- "IBD"
mds$QC_removal[match(pi_hat$IID2, mds$IID)] <- "IBD"

mds$QC_removal <- as.factor(mds$QC_removal)

plot(mds$C1, mds$C2, col=mds$QC_removal, xlab="First dimension", 
     ylab="Second dimension", pch=16)
plot(mds$C1, mds$C3, col=mds$QC_removal, xlab="First dimension", 
     ylab="Third dimension", pch=16)
plot(mds$C2, mds$C3, col=mds$QC_removal, xlab="Second dimension", 
     ylab="Third dimension", pch=16)
plot(mds$C3, mds$C4, col=mds$QC_removal, xlab="Third dimension", 
     ylab="Fourth dimension", pch=16)

mds[identify(mds$C3, mds$C4), ]
      FID   IID SOL        C1         C2         C3         C4         C5          C6          C7          C8         C9
100 FG106 FG106   0 0.0383056 0.00541788 0.00675128 -0.0329981 -0.0799293 -0.00327396 -0.00171556 -0.00136093 0.00222958
441  S306  S306   0 0.0374802 0.00516969 0.00602578 -0.0332917 -0.0805914 -0.00386855 -0.00294023 -0.00142025 0.00182714
567  S439  S439   0 0.0374608 0.00517456 0.00527530 -0.0321140 -0.0777346 -0.00528650 -0.00304928 -0.00129398 0.00114706
           C10 QC_removal
100 0.00342965       None
441 0.00294521       None
567 0.00162406       None

plot(mds$C4, mds$C5, col=mds$QC_removal, xlab="Fourth dimension", 
     ylab="Fifth dimension", pch=16)

mds[identify(mds$C4, mds$C5), ]
# Same

plot(mds$C5, mds$C6, col=mds$QC_removal, xlab="Fifth dimension", 
     ylab="Sixth dimension", pch=16)

mds[identify(mds$C5, mds$C6), ]
      FID   IID SOL        C1         C2         C3          C4          C5          C6          C7          C8          C9
100 FG106 FG106   0 0.0383056 0.00541788 0.00675128 -0.03299810 -0.07992930 -0.00327396 -0.00171556 -0.00136093  0.00222958
198   S55   S55   0 0.0659108 0.01518190 0.00552351 -0.00216074  0.01409360  0.01515300 -0.04677820 -0.03868940  0.04268600
313  S175  S175   0 0.0678426 0.01572220 0.00522490 -0.00263548  0.01421150  0.01532370 -0.04673760 -0.03753400  0.04261440
350  S214  S214   0 0.0408649 0.01374170 0.00166909  0.00298289  0.00493699 -0.03183100  0.05350580 -0.04609460 -0.02613780
351  S215  S215   0 0.0420967 0.01361930 0.00172898  0.00353858  0.00550267 -0.03172370  0.05354940 -0.04589050 -0.02700990
441  S306  S306   0 0.0374802 0.00516969 0.00602578 -0.03329170 -0.08059140 -0.00386855 -0.00294023 -0.00142025  0.00182714
567  S439  S439   0 0.0374608 0.00517456 0.00527530 -0.03211400 -0.07773460 -0.00528650 -0.00304928 -0.00129398  0.00114706
           C10 QC_removal
100 0.00342965       None
198 0.04177560       None
313 0.04133100       None
350 0.03844960       None
351 0.03861580       None
441 0.00294521       None
567 0.00162406       None

plot(mds$C7, mds$C8, col=mds$QC_removal, xlab="Seventh dimension", 
     ylab="Eight dimension", pch=16)

mds[identify(mds$C7, mds$C8), ]
11  FG11 FG11   0 0.0103025 -0.02364380 -0.015648500  0.000970638 0.000616365  0.00335050 -0.0190665  0.000280708 -0.03310950
174  S29  S29   0 0.0553443  0.01377090  0.003314680 -0.002572340 0.007797610  0.00157821  0.0192995  0.068826400  0.00727036
175  S30  S30   0 0.0528121  0.01320690  0.003151090 -0.002843760 0.008354150  0.00128666  0.0183132  0.069331600  0.00744520
198  S55  S55   0 0.0659108  0.01518190  0.005523510 -0.002160740 0.014093600  0.01515300 -0.0467782 -0.038689400  0.04268600
313 S175 S175   0 0.0678426  0.01572220  0.005224900 -0.002635480 0.014211500  0.01532370 -0.0467376 -0.037534000  0.04261440
350 S214 S214   0 0.0408649  0.01374170  0.001669090  0.002982890 0.004936990 -0.03183100  0.0535058 -0.046094600 -0.02613780
351 S215 S215   0 0.0420967  0.01361930  0.001728980  0.003538580 0.005502670 -0.03172370  0.0535494 -0.045890500 -0.02700990
482 S350 S350   0 0.0132485  0.00954125  0.000206884  0.007802540 0.002847410 -0.00929190  0.0158929 -0.006526010  0.00917633
483 S351 S351   0 0.0402828  0.01260200  0.001670620  0.003951040 0.005655830 -0.00896619  0.0157545 -0.004553960  0.00836644
            C10 QC_removal
11   0.00821355       None
174  0.04919020       None
175  0.04926670       None
198  0.04177560       None
313  0.04133100       None
350  0.03844960       None
351  0.03861580       None
482 -0.02374440       None
483 -0.02802640       None

plot(mds$C9, mds$C10, col=mds$QC_removal, xlab="Ninth dimension", 
     ylab="Tenth dimension", pch=16)

mds[identify(mds$C9, mds$C10), ]

     FID  IID SOL        C1          C2           C3           C4          C5          C6         C7           C8          C9
11  FG11 FG11   0 0.0103025 -0.02364380 -0.015648500  0.000970638 0.000616365  0.00335050 -0.0190665  0.000280708 -0.03310950
174  S29  S29   0 0.0553443  0.01377090  0.003314680 -0.002572340 0.007797610  0.00157821  0.0192995  0.068826400  0.00727036
175  S30  S30   0 0.0528121  0.01320690  0.003151090 -0.002843760 0.008354150  0.00128666  0.0183132  0.069331600  0.00744520
198  S55  S55   0 0.0659108  0.01518190  0.005523510 -0.002160740 0.014093600  0.01515300 -0.0467782 -0.038689400  0.04268600
313 S175 S175   0 0.0678426  0.01572220  0.005224900 -0.002635480 0.014211500  0.01532370 -0.0467376 -0.037534000  0.04261440
350 S214 S214   0 0.0408649  0.01374170  0.001669090  0.002982890 0.004936990 -0.03183100  0.0535058 -0.046094600 -0.02613780
351 S215 S215   0 0.0420967  0.01361930  0.001728980  0.003538580 0.005502670 -0.03172370  0.0535494 -0.045890500 -0.02700990
482 S350 S350   0 0.0132485  0.00954125  0.000206884  0.007802540 0.002847410 -0.00929190  0.0158929 -0.006526010  0.00917633
483 S351 S351   0 0.0402828  0.01260200  0.001670620  0.003951040 0.005655830 -0.00896619  0.0157545 -0.004553960  0.00836644
            C10 QC_removal
11   0.00821355       None
174  0.04919020       None
175  0.04926670       None
198  0.04177560       None
313  0.04133100       None
350  0.03844960       None
351  0.03861580       None
482 -0.02374440       None
483 -0.02802640       None

pc.outliers <- c("FG106", "S306", "S439", "S55", "S175", "S214", "S215", "FG11", "S29", "S30", "S350", "S351")
write.table(pc.outliers, "fail_pc_outliers.txt", sep="\t", quote=F, row.names=F,
col.names=F)
write.table(outliers, "PC-outliers.txt", sep="\t", quote=F, row.names=F)
```

Weirdly these seem to fall into pairs that aren't identified in the pi_hat analysis...


Kate identifies outliers: pop_outliers.txt

Compare the lists

```{r}
ke <- read.delim("../Genotyping/pop_outliers.txt", header=F)
kb <- read.delim("../Genotyping/pop_outliers_Katie.txt")
table(ke$V1 %in% kb$x)
table(kb$x %in% ke$V1)
```

# Remove selected individuals
# Repeat MDS again:

```{bash}
awk '{print $1,$1}' fail_pc_outliers.txt > fail_pc_outliers2.txt
cat fail-qc-inds2.txt fail_pc_outliers2.txt | sort -k1 | uniq > fail-qc-inds-pc1.txt
#120

~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked --remove fail-qc-inds-pc1.txt --extract AD_736_multi_ethnic_chip_updated_for_sex_check.prune.in --genome --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_first_round_ind_removed

~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked --remove fail-qc-inds-pc1.txt --read-genome AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_first_round_ind_removed.genome --cluster --mds-plot 10 --extract AD_736_multi_ethnic_chip_updated_for_sex_check.prune.in --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_first_round_ind_removed

```{r}
mds <- read.delim("AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_first_round_ind_removed.mds", sep="")
plot(mds$C1, mds$C2, xlab="First dimension", ylab="Second dimension", pch=16)
plot(mds$C1, mds$C3, xlab="First dimension", ylab="Third dimension", pch=16)
plot(mds$C2, mds$C3, xlab="Second dimension", ylab="Third dimension", pch=16)
plot(mds$C3, mds$C4, xlab="Third dimension", ylab="Fourth dimension", pch=16)
plot(mds$C4, mds$C5, pch=16)
plot(mds$C5, mds$C6, pch=16)
plot(mds$C6, mds$C7, pch=16)

mds[identify(mds$C6, mds$C7), ]
     FID  IID SOL        C1         C2          C3          C4          C5        C6        C7         C8          C9        C10
375 S247 S247   0 0.0142171 0.00304479 -0.00332059 -0.00464190 -0.00210783 0.0106935 0.0833229 0.00667391 -0.00419210 0.00172525
377 S249 S249   0 0.0160624 0.00353386 -0.00252363 -0.00529157 -0.00146149 0.0107407 0.0839351 0.00692201 -0.00329976 0.00172572

plot(mds$C7, mds$C8, pch=16)
mds[identify(mds$C7, mds$C8), ]
     FID  IID SOL           C1          C2          C3          C4           C5          C6           C7          C8          C9
227  S90  S90   0 -0.000299178 -0.00937543 -0.00821488  0.00131465 -0.000858114 -0.00878775 -0.000202076 -0.03117840 -0.06570900
228  S91  S91   0  0.000791149 -0.00891783 -0.00798274  0.00147183 -0.001010940 -0.00774210 -0.000228182 -0.03126160 -0.06644590
375 S247 S247   0  0.014217100  0.00304479 -0.00332059 -0.00464190 -0.002107830  0.01069350  0.083322900  0.00667391 -0.00419210
377 S249 S249   0  0.016062400  0.00353386 -0.00252363 -0.00529157 -0.001461490  0.01074070  0.083935100  0.00692201 -0.00329976
            C10
227 -0.00540909
228 -0.00520061
375  0.00172525
377  0.00172572

plot(mds$C8, mds$C9, pch=16)
mds[identify(mds$C8, mds$C9), ]
          FID       IID SOL           C1          C2           C3           C4           C5           C6           C7          C8
148        S2        S2   0 -0.000977771 -0.00575480 -0.007267260  0.000539431 -0.000389850 -0.006194820  0.000809842 -0.01304670
227       S90       S90   0 -0.000299178 -0.00937543 -0.008214880  0.001314650 -0.000858114 -0.008787750 -0.000202076 -0.03117840
228       S91       S91   0  0.000791149 -0.00891783 -0.007982740  0.001471830 -0.001010940 -0.007742100 -0.000228182 -0.03126160
322 S191_rep1 S191_rep1   0  0.006000340  0.00786130 -0.000676751  0.008482520 -0.000365551 -0.000284687 -0.000295478 -0.01806560
376      S248      S248   0 -0.001311690 -0.00394043 -0.005878800 -0.000585156 -0.001327020 -0.005490280  0.000380697 -0.00771863
557      S441      S441   0  0.015224600  0.01060450 -0.000227126  0.006806760 -0.001211350 -0.000953538 -0.000527724 -0.02101440
            C9         C10
148 -0.0229049 -0.00619638
227 -0.0657090 -0.00540909
228 -0.0664459 -0.00520061
322  0.0119784  0.01747820
376 -0.0133755 -0.00625212
557  0.0113814  0.02008780

plot(mds$C9, mds$C10, pch=16)
mds[identify(mds$C9, mds$C10), ]
          FID       IID SOL           C1          C2           C3           C4           C5           C6           C7          C8
148        S2        S2   0 -0.000977771 -0.00575480 -0.007267260  0.000539431 -0.000389850 -0.006194820  0.000809842 -0.01304670
227       S90       S90   0 -0.000299178 -0.00937543 -0.008214880  0.001314650 -0.000858114 -0.008787750 -0.000202076 -0.03117840
228       S91       S91   0  0.000791149 -0.00891783 -0.007982740  0.001471830 -0.001010940 -0.007742100 -0.000228182 -0.03126160
322 S191_rep1 S191_rep1   0  0.006000340  0.00786130 -0.000676751  0.008482520 -0.000365551 -0.000284687 -0.000295478 -0.01806560
334      S203      S203   0  0.010853600 -0.00195800  0.002307430 -0.000718097  0.003696580 -0.002584970 -0.001662860  0.00290136
349      S221      S221   0  0.060937300  0.02055010  0.006138290 -0.006033250 -0.003457940 -0.002940930 -0.002371110  0.01117110
557      S441      S441   0  0.015224600  0.01060450 -0.000227126  0.006806760 -0.001211350 -0.000953538 -0.000527724 -0.02101440
561      S445      S445   0  0.033183000  0.01724020  0.003293790  0.004048830 -0.000741285 -0.000312908  0.001006980 -0.01059430
             C9         C10
148 -0.02290490 -0.00619638
227 -0.06570900 -0.00540909
228 -0.06644590 -0.00520061
322  0.01197840  0.01747820
334  0.00431773 -0.03338860
349 -0.00714378  0.03104770
557  0.01138140  0.02008780
561  0.00223865  0.01943490
247 S111 S111   0 0.0109761 -0.00244329 0.00121961 -0.000480678 0.00330976 -0.00297694 -0.00185175 0.00230925 0.00456915
           C10
247 -0.0327135

pc.outliers <- c("S247", "S249", "S90", "S91", "S2", "S191_rep", "S248", "S441", "S203", "S221", "S445", "S111")
write.table(pc.outliers, "fail_pc_outliers_round2.txt", sep="\t", quote=F, row.names=F,
col.names=F)

```
```{bash}
awk '{print $1,$1}' fail_pc_outliers_round2.txt > fail_pc_outliersround2.txt
cat fail-qc-inds.txt fail_pc_outliersround2.txt | sort -k1 | uniq > fail-qc-inds-2.txt
#132

~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked --remove fail-qc-inds-2.txt --extract AD_736_multi_ethnic_chip_updated_for_sex_check.prune.in --genome --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_second_round_ind_removed

# And again
~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked --remove fail-qc-inds-2.txt --read-genome AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_second_round_ind_removed.genome --cluster --mds-plot 10 --extract AD_736_multi_ethnic_chip_updated_for_sex_check.prune.in --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_second_round_ind_removed


```{r}
mds <- read.delim("AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_second_round_ind_removed.mds", sep="")
plot(mds$C1, mds$C2, xlab="First dimension", ylab="Second dimension", pch=16)
plot(mds$C1, mds$C3, xlab="First dimension", ylab="Third dimension", pch=16)
plot(mds$C2, mds$C3, xlab="Second dimension", ylab="Third dimension", pch=16)
plot(mds$C3, mds$C4, xlab="Third dimension", ylab="Fourth dimension", pch=16)
plot(mds$C4, mds$C5, pch=16)
plot(mds$C5, mds$C6, pch=16)
plot(mds$C6, mds$C7, pch=16)
plot(mds$C7, mds$C8, pch=16)
plot(mds$C8, mds$C9, pch=16)
mds[identify(mds$C8, mds$C9, mds$IID), ]
     FID  IID SOL        C1        C2           C3           C4           C5          C6           C7        C8        C9
188  S49  S49   0 0.0346959 0.0147479  0.004961310 -0.000263269 -0.001809780 -0.00477760  0.000693952 0.0133758 0.0462327
415 S297 S297   0 0.0153320 0.0129726 -0.000166176 -0.006802530 -0.000694208 -0.00133150 -0.001821100 0.0193209 0.0170671
594 S490 S490   0 0.0647860 0.0241273  0.008989750 -0.006930650 -0.000590895 -0.00216824  0.005072810 0.0150292 0.0622120
           C10
188 0.01198710
415 0.00970496
594 0.02603540

plot(mds$C9, mds$C10, pch=16)
mds$IID[identify(mds$C9, mds$C10, mds$IID)]
more.outliers <- c("FG19", "FG78", "FG86", "FG88", "FG93", "S24",  "S27",  "S49",
"S75",  "S79",  "S80",  "S156", "S159", "S195", "S297", "S301", "S309", "S310", "S382", "S490")
write.table(more.outliers, "fail_pc_outliers_round3.txt", sep="\t", quote=F, row.names=F,
col.names=F)
```

awk '{print $1,$1}' fail_pc_outliers_round3.txt > fail_pc_outliersround3.txt
cat fail-qc-inds-2.txt fail_pc_outliersround3.txt | sort -k1 | uniq > fail-qc-inds-3.txt
#152

~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked --remove fail-qc-inds-3.txt --extract AD_736_multi_ethnic_chip_updated_for_sex_check.prune.in --genome --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_third_round_ind_removed

# And again
~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked --remove fail-qc-inds-3.txt --read-genome AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_third_round_ind_removed.genome --cluster --mds-plot 10 --extract AD_736_multi_ethnic_chip_updated_for_sex_check.prune.in --out AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_third_round_ind_removed
```

```r{}
mds <- read.delim("AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked_third_round_ind_removed.mds", sep="")
plot(mds$C1, mds$C2, xlab="First dimension", ylab="Second dimension", pch=16)
plot(mds$C1, mds$C3, xlab="First dimension", ylab="Third dimension", pch=16)
plot(mds$C2, mds$C3, xlab="Second dimension", ylab="Third dimension", pch=16)
plot(mds$C3, mds$C4, xlab="Third dimension", ylab="Fourth dimension", pch=16)
plot(mds$C4, mds$C5, pch=16)
plot(mds$C5, mds$C6, pch=16)
plot(mds$C6, mds$C7, pch=16)
plot(mds$C7, mds$C8, pch=16)
plot(mds$C8, mds$C9, pch=16)
plot(mds$C9, mds$C10, pch=16)
```

~/plink --bfile AD_736_multi_ethnic_chip_update_strand_update_parents_FID_sex_checked --remove fail-qc-inds-3.txt --make-bed --out AD_736_multi_ethnic_chip_updated_inds_removed

#################################################################################
# Per-marker QC

~/plink --bfile AD_736_multi_ethnic_chip_updated_inds_removed --missing --out AD_736_multi_ethnic_chip_updated_inds_removed

R
pdf("missing.hist.inds.removed.pdf")
missing <- read.table("AD_736_multi_ethnic_chip_updated_inds_removed.lmiss", header=T)
hist(missing[, 5], breaks=100, main="Histogram of SNP missingness")
dev.off()

~/plink --bfile AD_736_multi_ethnic_chip_updated_inds_removed --hardy --out AD_736_multi_ethnic_chip_updated_inds_removed

awk '{ if($9 < 0.00001) {print }}' AD_736_multi_ethnic_chip_updated_inds_removed.hwe | wc -l

awk '{ if($9 < 0.00001) {print }}' AD_736_multi_ethnic_chip_updated_inds_removed.hwe | awk '{ if($7 > $8) {print }}' | wc -l
awk '{ if($9 < 0.00001) {print }}' AD_736_multi_ethnic_chip_updated_inds_removed.hwe | awk '{ if($7 > $8) {print $2}}' > HWE_SNPs.txt

# 10028 SNPs were p-value < 0.00001.  Of these, 69 were due to excess heterozygosity and the rest were due to excess homozygosity.  It is unclear whether this is a reflection of consanguinity in the Abu Dhabi population without comparing to another population – eg European. Remove excess heterozygosity only (HWE_SNPs.txt)

~/plink --bfile AD_736_multi_ethnic_chip_updated_inds_removed --geno 0.05 --exclude HWE_SNPs.txt --maf 0.01 --make-bed --out AD_736_multi_ethnic_chip_updated_inds_snps_removed

64410 MB RAM detected; reserving 32205 MB for main workspace.

1779696 variants loaded from .bim file.
586 people (227 males, 359 females) loaded from .fam.
586 phenotype values loaded from .fam.
--exclude: 1779658 variants remaining.
Warning: At least 31 duplicate IDs in --exclude file.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 533 founders and 53 nonfounders present.
Calculating allele frequencies... done.
Warning: 256192 het. haploid genotypes present (see
AD_736_multi_ethnic_chip_updated_inds_snps_removed.hh ); many commands treat
these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.998731.
5122 variants removed due to missing genotype data (--geno).
876588 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
897948 variants and 586 people pass filters and QC.
Among remaining phenotypes, 114 are cases and 472 are controls.
--make-bed to AD_736_multi_ethnic_chip_updated_inds_snps_removed.bed +
AD_736_multi_ethnic_chip_updated_inds_snps_removed.bim +
AD_736_multi_ethnic_chip_updated_inds_snps_removed.fam ... done.


~/plink --bfile AD_736_multi_ethnic_chip_updated_inds_snps_removed --indep-pairwise 50 5 0.2 --out AD_736_multi_ethnic_chip_updated_inds_snps_removed_noLD

# Use the file “high-LD-regions.txt” from Anderson et al
# This is hg 17 so use liftover to hg19

~/plink --bfile AD_736_multi_ethnic_chip_updated_inds_snps_removed --exclude /well/jknight/Sepsis/Genotyping/high_LD_regions_hg19.txt --range --make-bed --out AD_736_multi_ethnic_chip_updated_inds_snps_removed_noLongRangeLD

~/plink --bfile AD_736_multi_ethnic_chip_updated_inds_snps_removed_noLongRangeLD --extract AD_736_multi_ethnic_chip_updated_inds_snps_removed_noLD.prune.in --pca --out AD_736_multi_ethnic_chip_updated_inds_snps_removed_noLD_noLD_genotyping_pca_clean

sh update_build.sh AD_736_multi_ethnic_chip_updated_inds_snps_removed Multi-EthnicGlobal_A1-b38.strand AD_736_multi_ethnic_chip_clean_b38


# Check PCs
```{r}
geno.pc <- read.table("AD_736_multi_ethnic_chip_updated_inds_snps_removed_noLD_noLD_genotyping_pca_clean.eigenvec", sep="", header=F, row.names=2)
geno.pc$V1 <- NULL
plot(geno.pc$V3, geno.pc$V4, pch=16)
plot(geno.pc$V5, geno.pc$V6, pch=16)
plot(geno.pc$V7, geno.pc$V8, pch=16)
rownames(geno.pc)[identify(geno.pc$V7, geno.pc$V8)]
plot(geno.pc$V9, geno.pc$V10, pch=16)
rownames(geno.pc)[identify(geno.pc$V9, geno.pc$V10)]
plot(geno.pc$V11, geno.pc$V12, pch=16)
rownames(geno.pc)[identify(geno.pc$V11, geno.pc$V12)]
# "S78"  "S369"
# "S181" "S238"
```
##################################################################################
# For eQTL

~/plink --bfile AD_736_multi_ethnic_chip_updated_inds_removed --keep Genotyping_IDs_for_eQTL.txt --geno 0.05 --exclude HWE_SNPs.txt --maf 0.01 --make-bed --out AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed

Options in effect:
  --bfile AD_736_multi_ethnic_chip_updated_inds_removed
  --exclude HWE_SNPs.txt
  --geno 0.05
  --keep Genotyping_IDs_for_eQTL.txt
  --maf 0.01
  --make-bed
  --out AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed

64410 MB RAM detected; reserving 32205 MB for main workspace.
1779696 variants loaded from .bim file.
586 people (227 males, 359 females) loaded from .fam.
586 phenotype values loaded from .fam.
--exclude: 1779658 variants remaining.
Warning: At least 31 duplicate IDs in --exclude file.
--keep: 470 people remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 442 founders and 28 nonfounders present.
Calculating allele frequencies... done.
Warning: 203568 het. haploid genotypes present (see
AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed.hh ); many commands
treat these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate in remaining samples is 0.998755.
4868 variants removed due to missing genotype data (--geno).
875050 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
899740 variants and 470 people pass filters and QC.
Among remaining phenotypes, 95 are cases and 375 are controls.
--make-bed to AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed.bed +
AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed.bim +
AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed.fam ... done.


~/plink --bfile AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed --indep-pairwise 50 5 0.2 --out AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed_noLD

# Use the file “high-LD-regions.txt” from Anderson et al
# This is hg 17 so use liftover to hg19

~/plink --bfile AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed --exclude /well/jknight/Sepsis/Genotyping/high_LD_regions_hg19.txt --range --make-bed --out AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed_noLongRangeLD

~/plink --bfile AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed_noLongRangeLD --extract AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed_noLD.prune.in --pca --out AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed_noLD_noLD_genotyping_pca_clean

sh update_build.sh AD_736_multi_ethnic_chip_updated_eQTL_inds_snps_removed Multi-EthnicGlobal_A1-b38.strand AD_736_multi_ethnic_chip_eQTL_genotyping_b38

~/plink --bfile AD_736_multi_ethnic_chip_eQTL_genotyping_b38 --recode A --out AD_736_multi_ethnic_chip_eQTL_genotyping_b38