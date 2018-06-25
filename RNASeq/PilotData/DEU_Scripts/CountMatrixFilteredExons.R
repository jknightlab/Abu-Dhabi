# Make count matrix for exon level data

# Count data have been generated on a individual sample level
# Take all htseq-count results and melt them in to one big dataframe

# Make a list of the full paths to the final feature counts for all samples (208)
basedir <- "/well/jknight/AbuDhabiRNA/Katie/DiffExonUsage/PilotData/"
s.info <- read.delim("/well/jknight/kate/RNAseq/P160233-REX/README.txt")
sample.names <- s.info$DataName
myfiles <- paste(basedir, sample.names, "_exons.txt", sep="")

# create a list to store single sample tables
DT <- list()

# read each file as array element of DT, keep gene id and count, and rename
# with sample name only
for (i in 1:length(myfiles) ) {
  DT[[myfiles[i]]] <- read.table(myfiles[i], header = F, stringsAsFactors = FALSE)
  cnts <- gsub("/well/jknight/AbuDhabiRNA/Katie/DiffExonUsage/PilotData/", "", myfiles[i])
  cnts <- gsub("_exons.txt", "", cnts)
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

####################################

# take all data rows to a new table
data.all <- data[grep("^ENS", rownames(data), perl=TRUE, invert=FALSE), ]

# final merged table
head(data.all, 3)

# write data to file
write.table(data.all, file = "/well/jknight/AbuDhabiRNA/Katie/DiffExonUsage/Pilot-raw-count-data-exons.txt", sep="\t")

# cleanup intermediate objects
rm(y, z, i, DT)

###################################

# filter to keep only exons with at least 3 reads in at least 9 samples
data.a <- as.matrix(data.all)
b <- which(data.a < 3, arr.ind=T)
data.a[b] <- 0
c <- which(data.a >= 3, arr.ind=T)
data.a[c] <- 1
d <- apply(data.a, 1, sum)
x.det <- which(d >= 9)
length(x.det)
data.filt <- data.all[x.det, ]
dim(data.filt)

write.table(data.filt, file = "/well/jknight/AbuDhabiRNA/Katie/DiffExonUsage/Pilot-filtered-count-data-exons.txt", sep="\t")