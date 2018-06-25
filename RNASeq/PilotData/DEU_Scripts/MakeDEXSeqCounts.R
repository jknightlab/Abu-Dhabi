# Make the exon count DEXSeqDataSet

exon.counts <- read.delim("/well/jknight/AbuDhabiRNA/Katie/DiffExonUsage/Pilot-filtered-count-data-exons.txt")
flattenedFile <- list.files("/well/jknight/AbuDhabiRNA/Katie/DiffExonUsage", 
                            pattern="gff$", full.names=TRUE)
basename(flattenedFile)
flattenedFile <- flattenedFile[2]
flattenedFile <- read.delim(flattenedFile, header=F)

flattenedFile <- flattenedFile[x.det, ]
gene.id <- as.character(flattenedFile$V9)
gene.id <- strsplit(gene.id, "*gene_id ")
gene.id <- unlist(gene.id)
gene.id <- gene.id[seq(2, length(gene.id), 2)]

sample.names <- colnames(exon.counts)

s.info <- read.delim("/well/jknight/kate/RNAseq/P160233-REX/README.txt")
s.info <- s.info[match(sample.names, s.info$DataName), ]

full.s.info <- read.delim("~/Abu-Dhabi/RNASeq/PilotData/PilotSampleInfoExpanded.txt",
                          stringsAsFactors=F, row.names=1)
full.s.info$PatientNumber <- as.character(full.s.info$PatientNumber)
full.s.info$Disease[full.s.info$PatientNumber == 115] <- "T2D"
sampleTable <- cbind(s.info, full.s.info[match(s.info$DataName, full.s.info$SampleId), ])

suppressPackageStartupMessages( library( "DEXSeq" ) )
dxd <- DEXSeqDataSet(exon.counts, sampleData=sampleTable, featureID=rownames(exon.counts),
                     groupID=gene.id, design= ~ SampleId + exon + SampleType:exon)
dxd
save(dxd, file="~/Abu-Dhabi/RNASeq/PilotData/ExonCountsFiltered.RData")