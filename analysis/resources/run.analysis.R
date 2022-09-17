#**********************Run Analysis on a cohort or all data*********************
run.analysis<-function(map,otu,taxaLevel,cohort,is.taxonomy){


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
  t.test.file.name<-paste0("output/",taxaLevel,"/",taxaLevel,"_t-test_",cohort,".txt")
  write.table(combined.tables,t.test.file.name,sep = "\t",quote = FALSE)

  adjustedp.t.test.names<-c("adj.p.vals.HC.T1","adj.p.vals.HC.T2","adj.p.vals.T1.T2")
  plots<-lapply(adjustedp.t.test.names,function(x) get.box.plots(map,otu,combined.tables,FDR=0.1,order.by.column =x,is.taxonomy))
  names(plots)<-adjustedp.t.test.names
  box.plot.file.name<-paste0(paste0("output/",taxaLevel,"/",taxaLevel,"_t-test_boxplot_ordered_",comparisons,"_",cohort,".pdf"))
  invisible(mapply(save.plots,plots,file.name=box.plot.file.name))

}


#**************Run Analysis for both cohorts (taxonomies or pathways)**********
run.analysis.all<-function(t,otu.path,meta.path,is.taxonomy=TRUE){

  taxaLevel<-capitalize(t)
  print(taxaLevel)

  #Load count table and metadata
  if(is.taxonomy){
    file.path<-sprintf(otu.path,t)
    file.norm<-load.data(meta.path,file.path,normalize = TRUE)
  } else{
    file.norm<-load.pathways(meta.path,otu.path)
  }

  map<-file.norm[[1]]
  otu<-file.norm[[2]]

  #Create output directory
  dir.name <- paste0('output/',taxaLevel)
  if(!dir.exists(dir.name)){
    dir.create(dir.name,recursive = TRUE,showWarnings = FALSE)
  }

  run.analysis(map,otu,taxaLevel,cohort = "All",is.taxonomy)
  run.analysis(map,otu,taxaLevel,cohort = "UNC",is.taxonomy)
  run.analysis(map,otu,taxaLevel,cohort = "Denver",is.taxonomy)
}


