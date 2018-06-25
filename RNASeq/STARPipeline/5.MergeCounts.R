#!/usr/bin/env Rscript

# Make count matrix

# Make a list of the full paths to the final feature counts for all samples (208)
basedir <- "Counts/"
s.info <- read.delim("mapping.info.txt", sep="", header=FALSE)
sample.names <- s.info$V1
myfiles <- paste(basedir, sample.names, ".counts.txt", sep="")

# create a list to store single sample tables
DT <- list()

# read each file as array element of DT, keep gene id and count, and rename
# with sample name only
for (i in 1:length(myfiles) ) {
  tryCatch({
    DT[[myfiles[i]]] <- read.table(myfiles[i], header = F, stringsAsFactors = FALSE)
    cnts <- gsub("Counts/", "", myfiles[i])
    cnts <- gsub(".counts.txt", "", cnts)
    colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
  }, error=function(e) print(myfiles[i]))
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


# inspect and look at the top row names
head(data)
tail(data)

# take all data rows to a new table
data.summary <- data[grep("__.*", rownames(data), perl=TRUE, invert=FALSE), ]
data.all <- data[-grep("__.*", rownames(data), perl=TRUE, invert=FALSE), ]

# final merged table
head(data.all, 3)

# write data to file
write.table(data.all, file = "Full_count_data.txt", sep="\t")

# cleanup intermediate objects
rm(y, z, i, DT)
