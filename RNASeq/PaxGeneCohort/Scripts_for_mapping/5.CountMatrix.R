# Make count matrix

# Make a list of the full paths to the final feature counts for all samples (208)
basedir <- "/well/jknight/AbuDhabiRNA/Katie/mapping/4.gene.counts/"
s.info <- read.delim("/well/jknight/AbuDhabiRNA/Katie/mapping/2.mapping/mapping.info.txt", 
                     sep="", header=FALSE)
sample.names <- s.info$V1
myfiles <- paste(basedir, sample.names, ".counts.txt", sep="")

# create a list to store single sample tables
DT <- list()

# read each file as array element of DT, keep gene id and count, and rename
# with sample name only
for (i in 1:length(myfiles) ) {
  DT[[myfiles[i]]] <- read.table(myfiles[i], header = T, stringsAsFactors = FALSE)
  DT[[myfiles[i]]] <- DT[[myfiles[i]]][, c(1, 7)]
  cnts <- gsub("/well/jknight/AbuDhabiRNA/Katie/mapping/4.gene.counts/", "", myfiles[i])
  cnts <- gsub(".counts.txt", "", cnts)
  colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
}

# merge all elements based on first ID columns
data <- DT[[myfiles[1]]]

# inspect
head(data)

# now add each additional table with the ID column as key
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

# write summary to file
write.table(data.all.summary, file = "~/Abu Dhabi/PilotRNA-hisat_counts_all-summary.txt", sep="\t")

####################################
# take all data rows to a new table
data.all <- data[grep("^ENS", rownames(data), perl=TRUE, invert=FALSE), ]

# final merged table
head(data.all, 3)

# write data to file
write.table(data.all, file = "~/Abu Dhabi/PaxGeneCohort/Full-count-data.txt", sep="\t")
write.table(data.all, file = "/well/jknight/AbuDhabiRNA/PaxGene-count-data.txt", sep="\t")

# cleanup intermediate objects
rm(y, z, i, DT)
