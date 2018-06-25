# Compare STAR to Tophat2 data

library(ggplot2)
library(reshape2)

counts <- read.delim("U:/Abu-Dhabi/RNASeq/PilotData/PilotRNA-hisat-count-data.txt")
star <- read.delim("U:/Abu-Dhabi/RNASeq/STARPipeline/STAR_full_count_data.txt")

counts <- counts[rownames(counts) %in% rownames(star), ]
star <- star[rownames(star) %in% rownames(counts), ]

counts <- counts[match(rownames(star), rownames(counts)), 
                 match(colnames(star), colnames(counts))]
all(rownames(counts) == rownames(star))
all(colnames(counts) == colnames(star))

count.m <- melt(counts)
star.m <- melt(star)

df.m <- cbind(count.m, star.m$value)
colnames(df.m) <- c("Sample", "Hisat2", "STAR")

pdf("CompareMapping.pdf", onefile=T)

for(i in 1:208){
  df <- subset(df.m, Sample == colnames(counts[i]))
  
  gg <- ggplot(df, aes(Hisat2, STAR, group=Sample)) +
    geom_point() +
    theme_bw() +
    geom_abline(intercept=0, slope=1) +
    ggtitle(label=df$Sample[1])
  print(gg)
}
dev.off()

totals <- cbind("Hisat2"=colSums(counts),
                "STAR"=colSums(star))
totals <- data.frame(totals)
ggplot(totals, aes(Hisat2, STAR)) +
  geom_point() +
  theme_bw() +
  geom_abline(intercept=0, slope=1)
