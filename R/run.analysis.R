require(Hmisc)

source("R/load.data.R")
source("R/sample.selection.R")
source("R/plot.R")
source("R/statistics.R")
source("R/save.plot.R")


#**********************Run Analysis on a cohort or all data*********************
run.analysis<-function(map,otu,taxaLevel,cohort){


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

# Principal Coordinate Analysis --------------------------------------------

  if(cohort!="All"){
    title=paste0(taxaLevel,"_",cohort)
  }else{
    title<-taxaLevel
  }
  p1<-pco.plot(map,otu,color = "Type", shape="Location" ,plot.title = title,show.legend = FALSE)
  pco.file.name<-paste0("output/",taxaLevel,"_pco_",cohort,".pdf")
  pdf(pco.file.name,width = 5,height = 5)
  print(p1)
  dev.off()

# Compare Groups ----------------------------------------------------------

  comparisons<-c("HC-T1","HC-T2","T1-T2")
  tables<-lapply(comparisons,function(x) compare.Groups(map,otu,group1 = strsplit(x,"-")[[1]][1],
                                                        group2 = strsplit(x,"-")[[1]][2],
                                                        paired = ifelse(strsplit(x,"-")[[1]][1]=="T1",TRUE,FALSE)))
  combined.tables<-cbind(tables[[1]],tables[[2]],tables[[3]])
  combined.tables<-combined.tables[,!duplicated(colnames(combined.tables))]
  t.test.file.name<-paste0("output/",taxaLevel,"_t-test_",cohort,".txt")
  write.table(combined.tables,t.test.file.name,sep = "\t",quote = FALSE)

  adjustedp.t.test.names<-c("adj.p.vals.HC.T1","adj.p.vals.HC.T2","adj.p.vals.T1.T2")
  plots<-lapply(adjustedp.t.test.names,function(x) get.box.plots(map,otu,combined.tables,order.by.column =x))
  names(plots)<-adjustedp.t.test.names
  box.plot.file.name<-paste0(paste0("output/",taxaLevel,"_t-test_boxplot_ordered_",comparisons,"_",cohort,".pdf"))
  mapply(save.plots,plots,file.name=box.plot.file.name)

# Mixed Linear Regression Analyses ----------------------------------------

  map.AN<-map %>% filter(Type != "HC")
  otu.AN<-otu %>% filter(map$Type != "HC")

  variables=c("BMI","BMIchange","BMIchangePerDay","TimeBetweenRecovery")
  variable.names<-c("BMI (kg/m2)","BMI Change (kg/m2)","BMI change per Day (kg/m2)","Time Between Recovery (d)")

  MLM.file.name<-paste0("output/",taxaLevel,"_MLM_",variables,"_",cohort,".txt")
  MLM.file.name.pdf<-paste0("output/",taxaLevel,"_MLM_",variables,"_",cohort,".pdf")
  MLM.result<-perform.MLM.all.vars(map.AN,otu.AN,variables,MLM.file.name)
  MLM.plots<-lapply(1:length(variables),function(x) get.scatter.plots(map.AN,otu.AN,variables[x],
                                                                      variable.names[x],result.test = MLM.result[[x]],legend.show = FALSE))
  mapply(save.plots,MLM.plots,file.name=MLM.file.name.pdf)

# Multivariate Analysis ---------------------------------------------------

  adonis.file.name<-paste0("output/",taxaLevel,"_Adonis_",cohort,".txt")
  adonis.file.name.pdf<-paste0("output/",taxaLevel,"_Adonis_",cohort,".pdf")
  adonis.result<-perform.adonis.all.vars(otu.AN,map.AN,variables,file.Name = adonis.file.name)
  adonis.plot<-plot.adonis(adonis.result,variable.names,adonis.file.name.pdf,show.legend = FALSE)
}


#**************Run Analysis for both cohorts (taxonomies or pathways)**********
run.analysis.all<-function(t,otu.path,meta.path,taxanomy=TRUE){
  taxaLevel<-capitalize(t)
  if(taxanomy){
    file.path<-sprintf(otu.path,t)
    file.norm<-load.data(meta.path,file.path,normalize = FALSE)
  } else{
    file.norm<-load.pathways(meta.path,otu.path)
  }

  map<-file.norm[[1]]
  otu<-file.norm[[2]]
  run.analysis(map,otu,taxaLevel,cohort = "All")
  run.analysis(map,otu,taxaLevel,cohort = "UNC")
  run.analysis(map,otu,taxaLevel,cohort = "Denver")
}


