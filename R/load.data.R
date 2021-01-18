load.data<-function(mapFile,otuFile,minPrevalence = 0.25,normalize=TRUE){

  map <- read.table(mapFile,sep="\t",header=TRUE,check.names = FALSE)
  otu<-read.table(otuFile,sep="\t",header=TRUE,check.names = FALSE,
                  row.names = 1, comment.char = "",quote = "")

  rownames(otu)<-sapply(rownames(otu),function(x){strsplit(x,"_")[[1]][1]})

  otu<-otu[map$SampleID,]

  if(normalize){

    averageReadPerSample = mean(rowSums(otu))
    otu_relab<-sweep(otu,1,rowSums(otu),"/")
    otu<-log10(otu_relab*averageReadPerSample + 1)
  }
    otu<-otu[,colMeans(otu>0)>minPrevalence]

    list(map,otu)

}


load.pathways<-function(map.file,pathway.file,minPrevalence = 0.25,unstratified=TRUE){
  map <- read.table(map.file,sep="\t",header=TRUE,check.names = FALSE)
  pathway<-read.table(pathway.file,sep="\t",header=TRUE,check.names = FALSE,
                  row.names = 1, comment.char = "",quote = "")

  if(unstratified)
    pathway<-pathway[-grep("|",rownames(pathway),fixed = TRUE),]

  pathway<-as.data.frame(t(pathway))
  rownames(pathway)<-sapply(rownames(pathway),function(x){strsplit(x,"_")[[1]][1]})
  pathway<-pathway[map$SampleID,]
  pathway<-pathway[,colMeans(pathway>0)>minPrevalence]

  list(map,pathway)
}


