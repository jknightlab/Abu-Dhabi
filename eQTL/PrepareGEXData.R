# Prepare gene expression data for eQTL mapping

library(limma)
library(edgeR)
library(DESeq2)
library(ggplot2)

# Raw count data for Paxgene samples
load("eQTL/dds_for_eQTL.Rdata")
dds
dds <- dds[mcols(dds)$chr %in% 1:22, ]
dds

load("U:/Abu-Dhabi/RNASeq/Full_transcript_metadata.RData")
colnames(mygtf)[1] <- "chr"
m1 <- match(rownames(assay(dds)), mygtf$gene_name)
assay.info <- mygtf[m1, ]
assay.info$TSS <- assay.info$start
assay.info$TSS[which(assay.info$strand == "-")] <- assay.info$end[which(assay.info$strand == "-")]

counts <- counts(dds)
s.info <- colData(dds)

dge <- DGEList(counts=counts, genes=assay.info, samples=s.info)

# Quantile normalisation and log2 transformation
counts.quantile <- voom(data.matrix(counts), plot=TRUE, normalize="quantile")

# PCA
pcaData <- prcomp(t(counts.quantile@.Data[[1]]), scale.=TRUE, center=TRUE)
ggplot(data.frame(pcaData$x), aes(x = PC1, y = PC2)) +
  geom_point(size=3) +
  coord_fixed() +
  theme_bw()

# Correct for total number of reads
counts.num <- apply(counts, 2, function(x) x/sum(x))

# PCA
pcaData <- prcomp(t(counts.num), scale.=TRUE, center=TRUE)
ggplot(data.frame(pcaData$x), aes(x = PC1, y = PC2)) +
  geom_point(size=3) +
  coord_fixed() +
  theme_bw()

# Remove outliers

# TMM and log2 transformation, correct for total mapped reads
dge.tmm <- calcNormFactors(dge, method="TMM")
# prior.count = average count to be added to each observation to avoid log(0)
dge.tmm.logcpm <- cpm(dge.tmm, log=TRUE, normalized.lib.sizes=TRUE, prior.count=1)

write.table(dge.tmm.logcpm, "eQTL/Normalised_GEX_data.txt", sep="\t", row.names=T)
write.table(assay.info, "U:/Abu-Dhabi/RNASeq/Full_transcript_metadata.txt", 
            sep="\t", row.names=F)
write.table(s.info, "U:/Abu-Dhabi/RNASeq/eQTL/sample_info.txt", sep="\t")
gex.final.t <- t(dge.tmm.logcpm)
write.table(gex.final.t, "eQTL/Normalised_GEX_data_t.txt", sep="\t", row.names=F)

# Generate gene expression PCs
gex.final.t <- read.delim("eQTL/Normalised_GEX_data_t.txt")
pcs <- prcomp(gex.final.t, scale=TRUE)
write.table(pcs$x, "eQTL/GEX_PCs.txt", sep="\t", quote=FALSE)

# Cumulative variance explained by the first 35 gene expression PCs
props <- (pcs$sdev^2)/sum(pcs$sdev^2) * 100

ggplot(data.frame(pcs$x), aes(x = PC1, y = PC2)) +
  geom_point(size=3) +
  coord_fixed() +
  theme_bw() +
  xlab(paste0("PC1 (", round(props[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(props[2], 2), "%)"))

ggplot(data.frame(pcs$x), aes(x = PC3, y = PC4)) +
  geom_point(size=3) +
  coord_fixed() +
  theme_bw() +
  xlab(paste0("PC3 (", round(props[3], 2), "%)")) +
  ylab(paste0("PC4 (", round(props[4], 2), "%)"))

for(x in seq(1, 19, 2)){
  p <- ggplot(data.frame(pcs$x), 
              aes_string(colnames(data.frame(pcs$x))[x], colnames(data.frame(pcs$x))[x + 1])) + 
    geom_point() + 
    theme_bw() + 
    xlab(paste("PC", x, "-", round(props[x], 2), "%")) + 
    ylab(paste("PC", x + 1, "-", round(props[x + 1], 2), "%"))  
  print(p)
  
}

props <- data.frame(props[1:50]*100)
colnames(props) <- "Props"
ggplot(props, aes(1:50, Props)) + geom_col() + theme_bw()

props <- cumsum(((pcs$sdev^2)/sum(pcs$sdev^2)))
props <- data.frame(props[1:100])
ggplot(props, aes(1:100, props.1.100.)) + geom_col() + theme_bw()
