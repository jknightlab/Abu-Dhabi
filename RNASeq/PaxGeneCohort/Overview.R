library(data.table)
library(DESeq2)
library(vsn)
library(ggplot2)
library(genefilter)
library(ggrepel)
library(reshape)
library(gridExtra)
library(gdata)

pdf("PaxGeneCohort/QC_plots.pdf", onefile=TRUE)

all.data <- read.delim("PaxGeneCohort/Full_count_data_STAR_filtered.txt")
to.rm <- which(colSums(all.data) == 0)
all.data <- all.data[, -to.rm]

mapping.stats <- read.delim("U:/Abu-Dhabi/RNASeq/PaxGeneCohort/STAR_mapping_stats.txt", 
                            row.names=NULL)
mapping.stats$Uniquelymappedreads. <- as.numeric(gsub("%", "", mapping.stats$Uniquelymappedreads.))
mapping.stats$X.ofreadsmappedtomultipleloci <- as.numeric(gsub("%", "", mapping.stats$X.ofreadsmappedtomultipleloci))
colnames(mapping.stats)[1] <- "Sample"

qc <- read.xls("U:/Abu-Dhabi/RNASeq/PaxGeneCohort/P170245_report_15Mar18.xlsx",
               sheet=3)
qc <- subset(qc, Sequenced.library. == "Yes")
mapping.stats <- mapping.stats[mapping.stats$Sample %in% qc$Sample.Name, ]

sample.info <- read.delim("PaxGeneCohort/CollatedFamilyInformation.txt",
                          stringsAsFactors=FALSE)

mapping.stats$Disease <- sample.info$Control.or.T2D[match(mapping.stats$Sample,
                                                          sample.info$Sample.name.on.plate)]
mapping.stats <- droplevels(mapping.stats)
counts <- all.data[, colnames(all.data) %in% qc$Sample.Name]

map.summary <- read.delim("U:/Abu-Dhabi/RNASeq/PaxGeneCohort/Summary_count_data_STAR_filtered.txt")
map.summary <- map.summary[, colnames(map.summary) %in% qc$Sample.Name]

cts <- data.frame(cbind("Mapped"=colSums(counts), "Unmapped"=colSums(map.summary)))

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

mapping.stats$Sample[which(mapping.stats$Uniquelymappedreads. < 60)]
# "S324" "S95"

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

k <- kOverA(k = round(ncol(counts) *0.05), A = 10)
ffun <- filterfun(k)
wh1 <- genefilter(counts, ffun)
sum(wh1) # check how many
filt <-  cbind(counts, wh1) # add column with TRUE/FALSE depending on pvalue filter
filt <- filt[ which(wh1 >= 1),] 
countDat <- filt[, -ncol(filt)]
dim(countDat) # 18682   620

sample.info <- sample.info[match(colnames(countDat), 
                                 sample.info$Sample.name.on.plate), ]
rownames(sample.info) <- colnames(countDat)
dds <- DESeqDataSetFromMatrix(countData=countDat,
                              colData=sample.info,
                              design= ~ Control.or.T2D)
dds

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(dds), 1, gm_mean)
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
vsd <- vst(dds, blind=TRUE)

plot(density(assay(vsd)[, 1]), ylim=c(0, 0.25), cex.main=0.8)
for (i in 2:ncol(assay(vsd))){
  lines(density(assay(vsd)[, i]))
}

d.tn <- dist(t(assay(vsd)))
plot(hclust(d.tn), cex=0.6, main="Cluster: VSN Data")
plot(hclust(d.tn), cex=0.6, main="Cluster: VSN Data", labels=vsd$RE)
# FG41F?

corMatrix <- cor(assay(vsd))
# write.table(corMatrix, "U:/Abu-Dhabi/RNASeq/PaxGeneCohort/corMatrix_norm.txt", sep="\t")

pcaData <- plotPCA(vsd, intgroup = c("Control.or.T2D", "RESEARCH.."), 
                   returnData = TRUE, ntop=18682)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, colour=Control.or.T2D)) +
  geom_point() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()

pcs <- prcomp(t(assay(vsd)), scale.=TRUE)
pcs <- data.frame(pcs$x)
pcs$Control.or.T2D <- vsd$Control.or.T2D

g1 <- ggplot(pcs, aes(x = PC1, y = PC2, colour=Control.or.T2D)) +
  geom_point() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()

g2 <- ggplot(pcs, aes(x = PC2, y = PC3, colour=Control.or.T2D)) +
  geom_point() +
  xlab(paste0("PC2: ", percentVar[2], "% variance")) +
  ylab(paste0("PC3: ", percentVar[3], "% variance")) +
  coord_fixed() +
  theme_bw()

g3 <- ggplot(pcs, aes(x = PC3, y = PC4, colour=Control.or.T2D)) +
  geom_point() +
  xlab(paste0("PC3: ", percentVar[3], "% variance")) +
  ylab(paste0("PC4: ", percentVar[4], "% variance")) +
  coord_fixed() +
  theme_bw()

g4 <- ggplot(pcs, aes(x = PC4, y = PC5, colour=Control.or.T2D)) +
  geom_point() +
  xlab(paste0("PC4: ", percentVar[4], "% variance")) +
  ylab(paste0("PC5: ", percentVar[5], "% variance")) +
  coord_fixed() +
  theme_bw()

grid.arrange(g1, g2, g3, g4, ncol=2)

pcaData <- pcaData[pcaData$RESEARCH.. %in% (names(table(dds$RESEARCH..)[table(dds$RESEARCH..) > 1])), ]
pcaData$RESEARCH.. <- as.character(pcaData$RESEARCH..)

ggplot(pcaData, aes(x = PC1, y = PC2, colour=Control.or.T2D)) +
  geom_text(aes(label=RESEARCH..)) +
  geom_line(aes(group=RESEARCH..)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()

pcaData <- plotPCA(vsd, intgroup = c("Control.or.T2D", "RESEARCH..", "Sex"), 
                   returnData = TRUE, ntop=18682/10)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, colour=Sex)) +
  geom_point() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()

load("U:/PANDIT/transcript_metadata.RData")
colnames(attributeData)[1] <- "chr"
m1 <- match(rownames(assay(vsd)), attributeData$gene_id)
assay.info <- attributeData[m1, ]
mcols(vsd) <- assay.info
table(mcols(vsd)$gene_biotype)
table(mcols(vsd)$chr)

vsd.nosex <- vsd[mcols(vsd)$chr != "X" & mcols(vsd)$chr != "Y", ]
pcaData <- plotPCA(vsd.nosex, intgroup = c("Control.or.T2D", "RESEARCH..", "Sex"), 
                   returnData = TRUE, ntop=18682/5)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, colour=Sex)) +
  geom_point() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()

vsd.sex <- vsd[mcols(vsd)$chr == "X" | mcols(vsd)$chr == "Y", ]
pcaData <- plotPCA(vsd.sex, intgroup = c("Control.or.T2D", "RESEARCH..", "Sex"), 
                   returnData = TRUE, ntop=18682/5)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, colour=Sex)) +
  geom_point() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()

pcaData$ToLabel <- FALSE
pcaData$ToLabel[pcaData$Sex == "F" & pcaData$PC1 > -5] <- TRUE
pcaData$ToLabel[pcaData$Sex == "M" & pcaData$PC1 < 15] <- TRUE
ggplot(pcaData, aes(Sex, PC1)) + 
  geom_boxplot() + 
  geom_point() +
  theme_bw() +
  geom_text_repel(data=subset(pcaData, ToLabel == TRUE), 
                  aes(label=name))

vsd.dc <- vsd[mcols(vsd)$chr != "X" & mcols(vsd)$chr != "Y", 
              vsd$Control.or.T2D == "Control" | vsd$Control.or.T2D == "T2D"]
pcaData <- plotPCA(vsd.dc, intgroup = c("Control.or.T2D", "RESEARCH..", "Sex"), 
                   returnData = TRUE, ntop=18682/10)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, colour=Control.or.T2D, shape=Sex)) +
  geom_point() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()

pcs <- prcomp(t(assay(vsd.dc)), scale.=TRUE)
pcs <- data.frame(pcs$x)
pcs$Control.or.T2D <- vsd.dc$Control.or.T2D

g1 <- ggplot(pcs, aes(x = PC1, y = PC2, colour=Control.or.T2D)) +
  geom_point() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()

g2 <- ggplot(pcs, aes(x = PC2, y = PC3, colour=Control.or.T2D)) +
  geom_point() +
  xlab(paste0("PC2: ", percentVar[2], "% variance")) +
  ylab(paste0("PC3: ", percentVar[3], "% variance")) +
  coord_fixed() +
  theme_bw()

g3 <- ggplot(pcs, aes(x = PC3, y = PC4, colour=Control.or.T2D)) +
  geom_point() +
  xlab(paste0("PC3: ", percentVar[3], "% variance")) +
  ylab(paste0("PC4: ", percentVar[4], "% variance")) +
  coord_fixed() +
  theme_bw()

g4 <- ggplot(pcs, aes(x = PC4, y = PC5, colour=Control.or.T2D)) +
  geom_point() +
  xlab(paste0("PC4: ", percentVar[4], "% variance")) +
  ylab(paste0("PC5: ", percentVar[5], "% variance")) +
  coord_fixed() +
  theme_bw()

grid.arrange(g1, g2, g3, g4, ncol=2)

save(dds, file="PaxGeneCohort/dds.Rdata")

dev.off()