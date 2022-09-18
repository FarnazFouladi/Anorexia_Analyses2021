set.seed(120)

source("../resources/load.packages.R")
source("../resources/load.data.R")
source("../resources/sample.selection.R")
source("../resources/plot.R")
source("../resources/statistics.R")
source("../resources/save.plot.R")
source("../resources/colors.R")


#**********************Run Analysis on a cohort or all data*********************
run.analysis<-function(map,otu,taxaLevel,cohort, outputSubDir, is.taxonomy){


  # Split Data --------------------------------------------------------------

  if(cohort=="UNC"){

    map.unc<-select.samples(map,otu,site = "UNC")[[1]]
    otu.unc<-select.samples(map,otu,site = "UNC")[[2]]
    map<-map.unc
    otu<-otu.unc


  } else if (cohort=="Denver"){

    map.denver<-select.samples(map,otu,site = "Denver")[[1]]
    otu.denver<-select.samples(map,otu,site = "Denver")[[2]]
    map<-map.denver
    otu<-otu.denver
  }

  if(cohort!="All"){
    title=paste0(taxaLevel,"_",cohort)
  }else{
    title<-taxaLevel
  }

  #Remove taxa with zero abundance
  otu <- otu[,colSums(otu) != 0]

  # Compare Groups ----------------------------------------------------------

  comparisons<-c("HC-T1","HC-T2","T1-T2")
  tables<-lapply(comparisons,function(x) compare.Groups(map,otu,group1 = strsplit(x,"-")[[1]][1],
                                                        group2 = strsplit(x,"-")[[1]][2],
                                                        paired = ifelse(strsplit(x,"-")[[1]][1]=="T1",TRUE,FALSE)))
  combined.tables<-cbind(tables[[1]],tables[[2]],tables[[3]])
  combined.tables<-combined.tables[,!duplicated(colnames(combined.tables))]
  t.test.file.name<-file.path(outputSubDir, paste0(taxaLevel,"_t-test_",cohort,".txt"))
  write.table(combined.tables,t.test.file.name,sep = "\t",quote = FALSE)

  adjustedp.t.test.names<-c("adj.p.vals.HC.T1","adj.p.vals.HC.T2","adj.p.vals.T1.T2")
  plots<-lapply(adjustedp.t.test.names,function(x) get.box.plots(map,otu,combined.tables,FDR=0.1,order.by.column =x,is.taxonomy))
  names(plots)<-adjustedp.t.test.names
  box.plot.file.name<-file.path(outputSubDir, paste0(taxaLevel,"_t-test_boxplot_ordered_",comparisons,"_",cohort,".pdf"))
  invisible(mapply(save.plots,plots,file.name=box.plot.file.name))
}


#**************Run Analysis for both cohorts (taxonomies or pathways)**********
run.analysis.all<-function(t,otu.path,meta.path, outputDir, is.taxonomy=TRUE){

  taxaLevel<-capitalize(t)
  print(taxaLevel)

  #Load count table and metadata
  if(is.taxonomy){
    file.path<-otu.path
    file.norm<-load.data(meta.path,file.path,normalize = TRUE)
  } else{
    file.norm<-load.pathways(meta.path,otu.path)
  }

  map<-file.norm[[1]]
  otu<-file.norm[[2]]

  #Create output subdirectory
  outputSubDir <- file.path(outputDir,taxaLevel)
  if(!dir.exists(outputSubDir)){
    dir.create(outputSubDir,recursive = TRUE,showWarnings = FALSE)
  }

  run.analysis(map,otu,taxaLevel,cohort = "All", outputSubDir, is.taxonomy)
  run.analysis(map,otu,taxaLevel,cohort = "UNC", outputSubDir, is.taxonomy)
  run.analysis(map,otu,taxaLevel,cohort = "Denver", outputSubDir, is.taxonomy)
}


