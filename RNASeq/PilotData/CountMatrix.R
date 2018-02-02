# Make count matrix

# Count data have been generated on a individual sample level
# Take all htseq-count results and melt them in to one big dataframe

# As input, the DESeq2 package expects count data as obtained, e.g., from RNA-seq or another high-throughput sequencing experiment, in the form of a matrix of integer values. The value in the i-th row and the j-th column of the matrix tells how many reads can be assigned to gene i in sample j. 

# The values in the matrix should be un-normalized counts or estimated counts of sequencing reads (for single-end RNA-seq) or fragments (for paired-end RNA-seq). It is important to provide count matrices as input for DESeq2’s statistical model (Love, Huber, and Anders 2014) to hold, as only the count values allow assessing the measurement precision correctly. The DESeq2 model internally corrects for library size, so transformed or normalized values such as counts scaled by library size should not be used as input.


# Make a list of the full paths to the final feature counts for all samples (208)
basedir <- "/well/jknight/kate/RNAseq/P160233-REX/"
s.info <- read.delim(paste(basedir, "README.txt", sep=""))
sample.names <- s.info$DataName
myfiles <- paste(basedir, sample.names, "_hisat/", "featureCount_dedup.txt", sep="")

# create a list to store single sample tables
DT <- list()

# read each file as array element of DT, keep gene id and count, and rename
# with sample name only
for (i in 1:length(myfiles) ) {
  DT[[myfiles[i]]] <- read.table(myfiles[i], header = T, stringsAsFactors = FALSE)
  DT[[myfiles[i]]] <- DT[[myfiles[i]]][, c(1, 7)]
  cnts <- gsub("/well/jknight/kate/RNAseq/P160233-REX/", "", myfiles[i])
  cnts <- gsub("_hisat/featureCount_dedup.txt", "", cnts)
  colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
}

# merge all elements based on first ID columns
data <- DT[[myfiles[1]]]

# inspect
head(data)

# ID ELL3153A180
# 1 ENSG00000223972           0
# 2 ENSG00000227232         232
# 3 ENSG00000243485           0
# 4 ENSG00000237613           1
# 5 ENSG00000268020           2
# 6 ENSG00000240361           0


# now add each other table with the ID column as key
for (i in 2:length(myfiles)) {
  y <- DT[[myfiles[i]]]
  z <- merge(data, y, by = c("ID"))
  data <- z
}

# ID column becomes rownames
rownames(data) <- data$ID
data <- data[, -1]

## add total counts per sample
data <- rbind(data, tot.counts=colSums(data))

# inspect and look at the top row names!
head(data)
tail(data)

####################################
# take summary rows to a new table
# ( not starting with ENS with invert=TRUE )
data.all.summary <- data[grep("^ENS", rownames(data), perl=TRUE, invert=TRUE), ]

# review
data.all.summary

# ELL3153A20 ELL3153A19 ELL3153A207 ELL3153A121 ELL3153A187
# ERCC-00002          0          0           0           0           0
# ERCC-00003          0          0           0           0           0
# ERCC-00004          0          0           0           0           0
# ERCC-00009          0          0           0           0           0
# ERCC-00012          0          0           0           0           0
# ERCC-00013          0          0           0           0           0
# ERCC-00014          0          0           0           0           0
# ERCC-00016          0          0           0           0           0
# ERCC-00017          0          0           0           0           0
# ERCC-00019          0          0           0           0           0
# ERCC-00022          0          0           0           0           0
# ERCC-00024          0          0           0           0           0
# ERCC-00025          0          0           0           0           0
# ERCC-00028          0          0           0           0           0
# ERCC-00031          0          0           0           0           0
# ERCC-00033          0          0           0           0           0
# ERCC-00034          0          0           0           0           0
# ERCC-00035          0          0           0           0           0
# ERCC-00039          0          0           0           0           0
# ERCC-00040          0          0           0           0           0
# ERCC-00041          0          0           0           0           0
# ERCC-00042          0          0           0           0           0
# ERCC-00043          0          0           0           0           0
# ERCC-00044          0          0           0           0           0
# ERCC-00046          0          0           0           0           0
# ERCC-00048          0          0           0           0           0
# ERCC-00051          0          0           0           0           0
# ERCC-00053          0          0           0           0           0
# ERCC-00054          0          0           0           0           0
# ERCC-00057          0          0           0           0           0
# ERCC-00058          0          0           0           0           0
# ERCC-00059          0          0           0           0           0
# ERCC-00060          0          0           0           0           0
# ERCC-00061          0          0           0           0           0
# ERCC-00062          0          0           0           0           0
# ERCC-00067          0          0           0           0           0
# ERCC-00069          0          0           0           0           0
# ERCC-00071          0          0           0           0           0
# ERCC-00073          0          0           0           0           0
# ERCC-00074          0          0           0           0           0
# ERCC-00075          0          0           0           0           0
# ERCC-00076          0          0           0           0           0
# ERCC-00077          0          0           0           0           0
# ERCC-00078          0          0           0           0           0
# ERCC-00079          0          0           0           0           0
# ERCC-00081          0          0           0           0           0
# ERCC-00083          0          0           0           0           0
# ERCC-00084          0          0           0           0           0
# ERCC-00085          0          0           0           0           0
# ERCC-00086          0          0           0           0           0
# ERCC-00092          0          0           0           0           0
# ERCC-00095          0          0           0           0           0
# ERCC-00096          0          0           0           0           0
# ERCC-00097          0          0           0           0           0
# ERCC-00098          0          0           0           0           0
# ERCC-00099          0          0           0           0           0
# ERCC-00104          0          0           0           0           0
# ERCC-00108          0          0           0           0           0
# ERCC-00109          0          0           0           0           0
# ERCC-00111          0          0           0           0           0
# ERCC-00112          0          0           0           0           0
# ERCC-00113          0          0           0           0           0
# ERCC-00116          0          0           0           0           0
# ERCC-00117          0          0           0           0           0
# ERCC-00120          0          0           0           0           0
# ERCC-00123          0          0           0           0           0
# ERCC-00126          0          0           0           0           0
# ERCC-00130          0          0           0           0           0
# ERCC-00131          0          0           0           0           0
# ERCC-00134          0          0           0           0           0
# ERCC-00136          0          0           0           0           0
# ERCC-00137          0          0           0           0           0
# ERCC-00138          0          0           0           0           0
# ERCC-00142          0          0           0           0           0
# ERCC-00143          0          0           0           0           0
# ERCC-00144          0          0           0           0           0
# ERCC-00145          0          0           0           0           0
# ERCC-00147          0          0           0           0           0
# ERCC-00148          0          0           0           0           0
# ERCC-00150          0          0           0           0           0
# ERCC-00154          0          0           0           0           0
# ERCC-00156          0          0           0           0           0
# ERCC-00157          0          0           0           0           0
# ERCC-00158          0          0           0           0           0
# ERCC-00160          0          0           0           0           0
# ERCC-00162          0          0           0           0           0
# ERCC-00163          0          0           0           0           0
# ERCC-00164          0          0           0           0           0
# ERCC-00165          0          0           0           0           0
# ERCC-00168          0          0           0           0           0
# ERCC-00170          0          0           0           0           0
# ERCC-00171          0          0           0           0           0
# gene0               0          0           0           0           0
# gene126             0          0           0           0           0
# gene132             0          0           0           0           0
# gene14              0          0           0           0           0
# gene144             0          0           0           0           0
# gene145             0          0           0           0           0
# gene15              0          0           0           0           0
# gene151             0          0           0           0           0
# gene17              0          0           0           0           0
# gene20              0          0           0           0           0
# gene22              0          0           0           0           0
# gene23              0          0           0           0           0
# gene24              0          0           0           0           0
# gene34              0          0           0           0           0
# gene4               0          0           0           0           0
# gene6               0          0           0           0           0
# gene67              0          0           0           0           0
# gene70              0          0           0           0           0
# gene71              0          0           0           0           0
# gene73              0          0           0           0           0
# rna12               0          0           0           0           0
# rna4                0          0           0           0           0
# rna5                0          0           0           0           0
# rna6                0          0           0           0           0
# rna9                0          0           0           0           0
# tot.counts    8831035   14629783    11860499    13972456    17237263

# write summary to file
write.table(data.all.summary, file = "~/Abu Dhabi/PilotRNA-hisat_counts_all-summary.txt", sep="\t")

####################################
# take all data rows to a new table
data.all <- data[grep("^ENS", rownames(data), perl=TRUE, invert=FALSE), ]

# final merged table
head(data.all, 3)

# write data to file
write.table(data.all, file = "~/Abu Dhabi/PilotRNA-hisat-count-data.txt", sep="\t")
write.table(data.all, file = "/well/jknight/AbuDhabiRNA/Pilot-count-data.txt", sep="\t")

# cleanup intermediate objects
rm(y, z, i, DT)
