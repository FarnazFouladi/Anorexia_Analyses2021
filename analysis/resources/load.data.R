#Load otu tables
load.data<-function(mapFile=meta.path,otuFile=file.path,minPrevalence = 0.25,normalize=TRUE){

  #Load meta data and otu table
  message("Loading meta data from file: ", mapFile)
  map <- read.table(mapFile,sep="\t",header=TRUE,check.names = FALSE,comment.char = "",quote = "")
  message("Loading OTU counts data from file: ", otuFile)
  otu<-read.table(otuFile,sep="\t",header=TRUE,check.names = FALSE,
                  row.names = 1, comment.char = "",quote = "")

  rownames(otu)<-sapply(rownames(otu),function(x){strsplit(x,"_")[[1]][1]})

  #Select samples that are found in the metadata
  otu<-otu[map$SampleID,]

  #Remove SK samples
  otu<-otu[map$Location!="SK",]
  map<-map[map$Location!="SK",]

  #Remove taxa with relative abundance less than 0.000001
  raThresh=0.000001
  message("Relative abundance threshold: ", raThresh)
  taxa_low_abun <- colSums(otu)/sum(colSums(otu)) < raThresh
  otu <- otu[,!taxa_low_abun]
  message("Removed ", sum(taxa_low_abun), " taxa.")

  if(normalize){
    averageReadPerSample = mean(rowSums(otu))
    otu_relab<-sweep(otu,1,rowSums(otu),"/")
    otu<-log10(otu_relab*averageReadPerSample + 1)
  }

  #Removing taxa that are present in less than 25% of samples.
  message("Removing taxa present is less than ", minPrevalence*100, "% of samples.")
  otu<-otu[,colMeans(otu>0)>minPrevalence]

  list(map,otu)
}

#Load pathway tables
load.pathways<-function(map.file,pathway.file,minPrevalence = 0.25,unstratified=TRUE){

  #Load metadata and metabolic pathways data
  message("Loading meta data from file: ", map.file)
  map <- read.table(map.file,sep="\t",header=TRUE,check.names = FALSE)
  message("Loading pathway data from file: ", pathway.file)
  pathway<-read.table(pathway.file,sep="\t",header=TRUE,check.names = FALSE,
                  row.names = 1, comment.char = "",quote = "")

  if(unstratified) pathway<-pathway[-grep("|",rownames(pathway),fixed = TRUE),]

  pathway<-as.data.frame(t(pathway))
  rownames(pathway)<-sapply(rownames(pathway),function(x){strsplit(x,"_")[[1]][1]})
  pathway<-pathway[map$SampleID,]

  #Remove SK samples
  pathway<-pathway[map$Location!="SK",]
  map<-map[map$Location!="SK",]

  #Remove pathways that are present in less than 25% of samples.
  message("Removing pathways present is less than ", minPrevalence*100, "% of samples.")
  pathway<-pathway[,colMeans(pathway>0)>minPrevalence]

  list(map,pathway)
}

load.ttest.results <- function(category, name){
  # Upstream modules
  taxaModuleName = dir(pipeRoot, pattern = "TaxonomyAnalysis", full.names = FALSE)
  pathwayModuleName = dir(pipeRoot, pattern = "PathwayAnalysis", full.names = FALSE)
  sourceModuleName = ifelse(category=="Pathway", pathwayModuleName, taxaModuleName)
  #
  file = file.path(pipeRoot, sourceModuleName, "output", category, paste0(category, name))
  message("Reading file: ", file)
  return(read.table(file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, quote = ""))
}

load.MLM.results <- function(category, name){
  # Upstream modules
  taxaModuleName = dir(pipeRoot, pattern = "TaxonomyAnalysis", full.names = FALSE)
  pathwayModuleName = dir(pipeRoot, pattern = "PathwayAnalysis", full.names = FALSE)
  sourceModuleName = ifelse(category=="Pathway", pathwayModuleName, taxaModuleName)
  #
  file = file.path(pipeRoot, sourceModuleName, "output", category, paste0(category, name))
  message("Reading file: ", file)
  return(read.table(file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, quote = ""))
}
