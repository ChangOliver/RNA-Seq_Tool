#!/usr/bin/RScript
# htseq-combine_all.R
# copy this code to RStudio and adapt file locations to match yours

# Take 'all' htseq-count results and melt them in to one big dataframe

# where are we?

print("Merging counts into one file...")
args <- commandArgs(trailingOnly = TRUE)
cntdir <- args[[1]]

pat <- "*.htseq-count.txt"
tophat.all <- list.files(path = cntdir,
                         pattern = pat,
                         all.files = TRUE,
                         recursive = FALSE,
                         ignore.case = FALSE,
                         include.dirs = FALSE)

# we choose the 'all' series
myfiles <- tophat.all
DT <- list()

# read each file as array element of DT and rename the last 2 cols
# we created a list of single sample tables
for (i in 1:length(myfiles) ) {
  infile = paste(cntdir, myfiles[i], sep = "/")
  DT[[myfiles[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
  cnts <- gsub("(.*).htseq-count.txt", "\\1", myfiles[i])
  colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
}

# merge all elements based on first ID columns
data <- DT[[myfiles[1]]]

# we now add each other table with the ID column as key
for (i in 2:length(myfiles)) {
  y <- DT[[myfiles[i]]]
  z <- merge(data, y, by = c("ID"))
  data <- z
}

# ID column becomes rownames
rownames(data) <- data$ID
data <- data[,-1]

## add total counts per sample
data <- rbind(data, tot.counts=colSums(data))

# take all data rows to a new table
data.all <- data[grep("^ENS", rownames(data), perl=TRUE, invert=FALSE), ]

# write data to file
write.csv(data.all, file = paste0(cntdir, "htseq_counts_all.csv"))

# cleanup intermediate objects
rm(y, z, i, DT)