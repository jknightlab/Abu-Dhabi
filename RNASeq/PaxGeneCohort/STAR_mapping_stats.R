# Mapping QC STAR AD Paxgene

library(ggplot2)
library(DESeq2)
library(reshape)
library(ggrepel)
library(gridExtra)
library(vsn)
library(gdata)

mapping.stats <- read.delim("U:/Abu-Dhabi/RNASeq/PaxGeneCohort/STAR_mapping_stats.txt", 
                            row.names=NULL)
mapping.stats$Uniquelymappedreads. <- as.numeric(gsub("%", "", mapping.stats$Uniquelymappedreads.))
mapping.stats$X.ofreadsmappedtomultipleloci <- as.numeric(gsub("%", "", mapping.stats$X.ofreadsmappedtomultipleloci))
colnames(mapping.stats)[1] <- "Sample"

qc <- read.xls("U:/Abu-Dhabi/RNASeq/PaxGeneCohort/P170245_report_15Mar18.xlsx",
               sheet=3)
qc <- subset(qc, Sequenced.library. == "Yes")
mapping.stats <- mapping.stats[mapping.stats$Sample %in% qc$Sample.Name, ]

sample.info <- read.xls("N:/jknight/Group_members_personal_data_storage/Lawrence/Abu Dabi Patient medical data/Collated Family information.xlsm")

mapping.stats$Disease <- sample.info$Control.or.T2D[match(mapping.stats$Sample,
                                                          sample.info$Sample.name.on.plate)]
mapping.stats <- droplevels(mapping.stats)

counts <- read.delim("U:/Abu-Dhabi/RNASeq/PaxGeneCohort/Full_count_data_STAR_filtered.txt")
counts <- counts[, colnames(counts) %in% qc$Sample.Name]
map.summary <- read.delim("U:/Abu-Dhabi/RNASeq/PaxGeneCohort/Summary_count_data_STAR_filtered.txt")
map.summary <- map.summary[, colnames(map.summary) %in% qc$Sample.Name]

cts <- data.frame(cbind("Mapped"=colSums(counts), "Unmapped"=colSums(map.summary)))

dds <- DESeqDataSetFromMatrix(countData=counts,
                               colData=mapping.stats,
                               design= ~ Disease)
dds
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]
# dds

pdf("U:/Abu-Dhabi/RNASeq/PaxGeneCohort/STARMappingQC.pdf", onefile=TRUE)

ggplot(mapping.stats, aes(Numberofinputreads)) +
  geom_histogram() +
  theme_bw()

ggplot(mapping.stats, aes(Uniquelymappedreadsnumber)) +
  geom_histogram() +
  theme_bw()

ggplot(mapping.stats, aes(Averagemappedlength)) +
  geom_histogram() +
  theme_bw()

ggplot(mapping.stats, aes(Numberofinputreads, Uniquelymappedreadsnumber, 
                          colour=X.ofreadsmappedtomultipleloci)) +
  scale_color_continuous(low="darkblue", high="red") +
  geom_point() +
  theme_bw() +
  xlab(label="Number of reads") +
  ylab(label="Number of uniquely mapping reads")

ggplot(mapping.stats, aes(Disease, Uniquelymappedreadsnumber)) +
  geom_boxplot() +
  geom_point() +
  theme_bw() +
  xlab(label="Sample Type") +
  ylab(label="Number of uniquely mapped reads")


ggplot(mapping.stats, aes(Numberofinputreads, Numberofsplices.Total)) +
  geom_point() +
  theme_bw() +
  xlab(label="Number of reads") +
  ylab(label="Number of splices identified")

ggplot(mapping.stats, aes(Numberofinputreads, Mismatchrateperbase..)) +
  geom_point() +
  theme_bw() +
  xlab(label="Number of reads") +
  ylab(label="Mismatch rate per base")

ggplot(mapping.stats, aes(Numberofinputreads, Numberofreadsmappedtomultipleloci)) +
  geom_point() +
  theme_bw() +
  xlab(label="Number of reads") +
  ylab(label="Number of multi-mapping reads")

ggplot(mapping.stats, aes(Numberofinputreads, Uniquelymappedreads.)) +
  geom_point() +
  theme_bw() +
  xlab(label="Number of reads") +
  ylab(label="Percentage uniquely mapped reads")

ggplot(mapping.stats, aes(Disease, Uniquelymappedreads.)) +
  geom_boxplot() +
  geom_point() +
  theme_bw() +
  xlab(label="Sample Type") +
  ylab(label="Percentage uniquely mapped reads")

ggplot(cts, aes(Mapped, Unmapped)) +
  geom_point() +
  theme_bw() +
  xlab(label="Total counted reads") +
  ylab(label="Total unmapped reads")

vsd <- vst(dds, blind=TRUE)

# Compare before and after
p1 <- meanSdPlot(assay(dds), plot=FALSE)
p1 <- p1$gg + theme_bw()
p2 <- meanSdPlot(assay(vsd), plot=FALSE)
p2 <- p2$gg + theme_bw()

grid.arrange(p1, p2, nrow=1)

melt.exprs <- melt(assay(vsd))
melt.exprs$Batch <- mapping.stats$Batch[match(melt.exprs$X2, 
                                              mapping.stats$Sample)]
ggplot(melt.exprs, aes(X2, value, colour=Batch)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position="bottom") +
  xlab("Sample")

pcaData <- plotPCA(vsd, intgroup = c("Batch"), returnData = TRUE,
                   ntop=35258)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()

dev.off()
