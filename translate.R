suppressMessages(library(httr))
suppressMessages(library(jsonlite))

get.translation <- function(res, col){
  idSymbolTable <- read.csv(file = "./hugo.txt", sep = '\t')
  hugoDict <- as.vector(idSymbolTable$Approved.symbol)
  names(hugoDict) <- idSymbolTable$Ensembl.gene.ID
  message(paste("Total gene:", nrow(res)))
  
  for (i in c(1:nrow(res))){
    if ((i+1) %% 5000 == 0){
      message(paste(i+1, "gene translated"))
    }
    
    id <- sub( "(ENS[A-Z]+[0-9]{11}).*", "\\1", res[i, 1])

    if (id %in% names(hugoDict)){
      res[i, col] <- hugoDict[[id]]
    }else{
      url = paste("https://biotools.fr/human/ensembl_symbol_converter/?api=1&id=", id, sep="")
      r <- GET(url)
      output = fromJSON(content(r, "text", encoding = "UTF-8"), flatten=TRUE)

      res[i, col] <- ifelse((length(output) != 0 & !is.null(output[[id]])), output[[id]], res[i, 1])
    }
  }
  
  return(res)
}

print("Translating Ensembl to Symbol...(This may take a while)")

args <- commandArgs(trailingOnly = TRUE)

cntDir <- args[[1]] #"../6.STATS/htseq_counts_all.csv"
counts <- paste0(cntDir, "htseq_counts_all.csv")
ctsTable <- read.table(counts, check.names=FALSE, header=TRUE, sep = ",", stringsAsFactors = FALSE)

df <- data.frame(as.matrix(ctsTable[,1]), stringsAsFactors = FALSE)
colnames(df)[1] <- "Ensembl.ID"
df[ , "Gene.name"] <- NA
df <- get.translation(df, 2)

print("Dumping translation table")
write.csv(df, "./translation_table.csv", row.names = FALSE)