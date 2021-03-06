#############################
# Load libraries
#############################
print("Loading libraries...")
suppressMessages(library(dplyr))
suppressMessages(library(DESeq2))
suppressMessages(library(OUTRIDER))
suppressMessages(library(IHW))
suppressMessages(library(tibble))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(rlist))
suppressMessages(library(magrittr))
suppressMessages(library(httr))
suppressMessages(library(jsonlite))

#############################
# Function declaration
#############################
get.translation <- function(res, col){
  idSymbolTable <- read.csv(file = paste0(srcDir, "/translation_table.csv"))
  translationDict <- as.vector(idSymbolTable$Gene.name)
  names(translationDict) <- idSymbolTable$Ensembl.ID
  
  for (i in c(1:nrow(res))){
    res[i, col] <- translationDict[[res[[i, 1]]]]
  }
  
  return(res)
}
get.count <- function(dir){
  
  #### Get files with suffix htseq-count.txt
  files <- grep("*.htseq-count.txt",list.files(dir),value=TRUE)
  
  #### Extract id
  names <- sub("RNA-(.*).htseq-count.txt", "\\1", files)
  
  #### Extract control index
  control <- grep("RNA-Control-.*", files, value=FALSE)
  
  #### Mark case or control
  condition <- rep("Case", times=length(files)) #names
  condition[control] <- "Control"
  
  #### create file info table, set conditions as factor
  table <- data.frame(sampleName = names, fileName = files, condition = condition)
  table$condition <- relevel(table$condition, ref = "Control")  
  return(list("table" = table, "condition" = condition, "names" = names))
}
run.DESeq2 <- function(dir, table, condition, resDir){

  #### Create DESeq object
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = table, directory = dir, design= ~ condition)
  
  #### Run differential expression and save dds object
  ddsHTSeq <- DESeq(ddsHTSeq, minReplicatesForReplace=Inf)
  saveRDS(ddsHTSeq, file=paste0(resDir, '/dds.rds')) 
  return(ddsHTSeq)
}
run.OUTRIDER <- function(countDir, resDir){
  #### Load data
  counts <- paste0(countDir, "/htseq_counts_all.csv")
  ctsTable <- read.table(counts, check.names=FALSE, header=TRUE, sep = ",", stringsAsFactors = FALSE)
  ctsMatrix <- as.matrix(ctsTable[,-1])
  rownames(ctsMatrix) <- as.character(get.translation(as.matrix(ctsTable[,1]), 1))
  
  #### Create OUTRIDER object
  ods <- OutriderDataSet(countData=ctsMatrix)
  ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
  ods <- OUTRIDER(ods)
  saveRDS(ods, file = paste0(resDir, "/ods.rds"))
  return(ods)
}
extract.result <- function(dds, mode, case){
  
  if (mode == "DESeq2"){
    #### Extract result
    res <- DESeq2::results(dds, contrast=c("condition", case, "Control"), cooksCutoff=FALSE, independentFiltering=FALSE)
  
    res <- data.frame(res, row.names = res@rownames)
    res <- subset(res, !is.na(res$padj))
    res[ , "name"] <- NA
    res <- rownames_to_column(res[c(7, 1:6)]) %>% get.translation(2)
  }
  else if (mode == "OUTRIDER"){
    res <- OUTRIDER::results(ods, all = TRUE)
    #res[ , "name"] <- as.character(NA)
    #res <- res[, c(1, 15, 2:14)] %>% get.translation(2)
    #colnames(res)[2] <- "geneID"
    #colnames(res)[1] <- "ensemblID"
  }
  
  return(res)
}

plot.volcano <- function(resList, resDir, case, fcCutoff, padjCutoff){
  
  imgPath <- paste0(getwd(), "/imgs/")
  dir.create(imgPath, showWarnings = FALSE)

  resList <- mutate(resList, 
                         group = ifelse( (log2FoldChange >= fcCutoff & -log10(padj) > padjCutoff), "R", 
                                            ifelse( (log2FoldChange <= (-fcCutoff) & -log10(padj) > padjCutoff), "G", "B")))
  
  resList <- resList[order(resList$group, decreasing = TRUE), ]
  volcano <- ggplot(resList, mapping = aes(x = log2FoldChange, y = -log10(padj))) + 
    geom_point(aes(color = group)) + 
    scale_color_manual(values = c("R" = "#EC7063", "G" = "#48C9B0", "B" = "#85929E")) + #1:red 2:green 3:gray
    ggtitle(case) +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") + 
    theme(legend.position = "none", 
          plot.title = element_text(size = rel(1.5), hjust = 0.5), 
          axis.title = element_text(size = rel(1.25)))
  ggsave(paste0(imgPath, "volcano.png"), width = 9.93, height = 6.86, dpi = 300)
    
  label_volcano <- volcano + geom_label_repel(data = filter(resList, group != "B"), 
                                              aes(label = as.character(name)), 
                                              size = 3) 
  ggsave(paste0(imgPath, "label_volcano.png"), width = 9.93, height = 6.86, dpi =300)
  
  return(resList)
}

plot_QQ_ExRank <- function(res, resDir){
  
  imgPath <- paste0(resDir, "/imgs/")
  dir.create(imgPath, showWarnings = FALSE)
  
  res <- subset(res, res$aberrant)
  
  for (i in c(1:nrow(res))){
    
      png(paste(paste0(imgPath, "QQ"), res[i]$sampleID, res[i]$geneID, ".png", sep = "_"), width = 715, height = 494)
      QQ <- plotQQ(ods, res[i, geneID])
      dev.off()
    
      # expression rank of a gene with outlier events
      EX <- plotExpressionRank(ods, res[i, geneID], basePlot=TRUE)
      ggsave(paste(paste0(imgPath, "EX"), res[i]$sampleID, res[i]$geneID, ".png", sep = "_"), width = 9.93, height = 6.86, dpi = 300)
  }
}


#############################
# Main
#############################
srcDir <- getwd()

print("Starting analysis")
args <- commandArgs(trailingOnly = TRUE)

mode <- "-all"#args[[1]]
dataPath <- "../6.STATS"#args[[2]]
resDir <- "../results/group"#args[[3]]
if (length(args) > 3){
  fcCutoff <- as.numeric(args[[4]])
  padjCutoff <- as.numeric(args[[5]])
}else{
  fcCutoff <- 2.5
  padjCutoff <- 2
}

dir.create(resDir, showWarnings = FALSE)

task <- vector()
if (mode == "-all"){
	task <- c("DESeq2", "OUTRIDER")
} else if (mode == "-outrider"){
	task <- c("OUTRIDER")
} else if (mode == "-deseq2"){
	task <- c("DESeq2")
}

if ("DESeq2" %in% task){

	### Run DESeq2
	if (file_test("-d", dataPath)){
		### Fetch htseq data
	  print("Running DESeq2 analysis...")
		returnVal <- get.count(dataPath)
		dds <- run.DESeq2(dataPath, returnVal$table, returnVal$condition, resDir)
	} else{
	  print("Reading rds data...")
		dds <- readRDS(file=dataPath)
	}
  
  setwd(resDir)
  resDir <- getwd()
  case_list <- seq(1,37)
  case_list <- sprintf("%03dA", case_list)
  for (i in case_list){
    print(paste0("Processing ", i))
    
    dir.create(paste0(resDir, "/RNA-", i), showWarnings = FALSE)
    setwd(paste0(resDir, "/RNA-", i))
    dir.create(paste0(getwd(), "/DESeq2"), showWarnings = FALSE)
    setwd(paste0(getwd(), "/DESeq2"))
    
  	message("extracting result")
  	desResult <- extract.result(dds, "DESeq2", i)
  
  	message("generating volcano plot")
  	desResult <- plot.volcano(desResult, resDir, i, fcCutoff, padjCutoff)
  
  	message("Dumping results...")
  	write.csv(desResult, paste0(getwd(), "/DESeq2_result.csv"), row.names = FALSE)
  	setwd(resDir)
  }
}

if ("OUTRIDER" %in% task){
	#### Run OUTRIDER
	print("Running OUTRIDER analysis...")
	if (file_test("-d", dataPath)){
		ods <- run.OUTRIDER(dataPath, resDir)
  } else{
		ods <- readRDS(file=dataPath)
  }

	message("extracting result")
	odsResult <- extract.result(ods, "OUTRIDER")

	message("generating QQ & ExRank plot")
	plot_QQ_ExRank(odsResult, resDir)

	print("Dumping results...")
	write.csv(odsResult, paste0(resDir, "/OUTRIDER_result.csv"), row.names = FALSE)
}

