#**********************Run Analysis on a cohort or all data*********************
run.analysis<-function(map,otu,taxaLevel,cohort,is.taxonomy){

variables=c("BMI")
variables.names<-c("BMI (kg/m2)")

variables.change=c("BMIchange","BMIchangePerDay")
variables.change.names<-c("BMI Change (kg/m2)","BMI change per Day (kg/m2)")


# Split Data --------------------------------------------------------------

  if(cohort=="UNC"){

    map.unc<-select.samples(map,otu,site = "UNC")[[1]]
    otu.unc<-select.samples(map,otu,site = "UNC")[[2]]
    map<-map.unc
    otu<-otu.unc

    variables<-c("BMI","TotalFatPercent","StandardBIA_FatPercent",
                 "Trunk..Percent.Fat","Head..Percent.Fat",
                 "Trunk..Fat.mass..g.","Head..Fat.mass..g.","Total..Fat.mass..g.")
    variables.names<-c("BMI (kg/m2)","Total Fat (percent)","Standard BIA Fat (percent)",
                      "Trunk Fat (percent)","Head Fat (percent)",
                      "Trunk Fat Mass (g)","Head Fat Mass (g)","Total Fat Mass (g)")
    variables.change<-c("BMIchange","BMIchangePerDay",
                 "trunkFatChange","headFatChange",
                 "totalFatChange","trunkPercentFatChange","headPercentFatChange",
                 "totalPercentFatChange")
    variables.change.names<-c("BMI Change (kg/m2)","BMI change per Day (kg/m2)",
                      "Trunk Fat Change (g)","Head Fat Change (g)", "Total Fat Change (g)",
                      "Trunk Fat Change (percent)","Head Fat Change (percent)", "Total Fat Change (percent)")


  } else if (cohort=="Denver"){

    map.denver<-select.samples(map,otu,site = "Denver")[[1]]
    otu.denver<-select.samples(map,otu,site = "Denver")[[2]]
    map<-map.denver
    otu<-otu.denver
    variables<-c("BMI")
    variables.names<-c("BMI (kg/m2)")
    variables.change<-c("BMIchange","BMIchangePerDay")
    variables.change.names<-c("BMI Change (kg/m2)","BMI change per Day (kg/m2)")
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
  mapply(save.plots,plots,file.name=box.plot.file.name)


# Baseline Microbiome -----------------------------------------------------
  map.AN<-map %>% filter(Type != "HC")
  otu.AN<-otu %>% filter(map$Type != "HC")
  otu.AN <- otu.AN[, colSums(otu.AN)!=0]

  map.AN.baseline <- map.AN %>% filter(Type == "T1")
  otu.AN.baseline <- otu.AN %>% filter(map.AN$Type == "T1")
  otu.AN.baseline <- otu.AN.baseline[, colSums(otu.AN.baseline)!=0]

# Mixed Linear Regression Analyses ----------------------------------------

  ## For variables including T1 and T2 microbiome

  MLM.file.name<-paste0("output/",taxaLevel,"/",taxaLevel,"_MLM_",variables,"_",cohort,".txt")
  MLM.file.name.pdf<-paste0("output/",taxaLevel,"/",taxaLevel,"_MLM_",variables,"_",cohort,".pdf")
  MLM.result<-perform.MLM.all.vars(map.AN,otu.AN,variables,MLM.file.name,changeInVariable = FALSE)
  MLM.plots<-lapply(1:length(variables),function(x) get.scatter.plots(map.AN,otu.AN,variables[x],
                                                                      variables.names[x],result.test = MLM.result[[x]],legend.show = FALSE))
  mapply(save.plots,MLM.plots,file.name=MLM.file.name.pdf)

  ## For change in variables including only T1 microbiome
  MLM.file.name<-paste0("output/",taxaLevel,"/",taxaLevel,"_MLM_",variables.change,"_",cohort,".txt")
  MLM.file.name.pdf<-paste0("output/",taxaLevel,"/",taxaLevel,"_MLM_",variables.change,"_",cohort,".pdf")
  MLM.result<-perform.MLM.all.vars(map.AN.baseline,otu.AN.baseline,variables.change,MLM.file.name,changeInVariable = TRUE)
  MLM.plots<-lapply(1:length(variables.change),function(x) get.scatter.plots(map.AN.baseline,otu.AN.baseline,variables.change[x],
                                                                      variables.change.names[x],result.test = MLM.result[[x]],legend.show = FALSE))
  mapply(save.plots,MLM.plots,file.name=MLM.file.name.pdf)


# Multivariate Analysis ---------------------------------------------------

  #For variables
  adonis.file.name<-paste0("output/",taxaLevel,"/",taxaLevel,"_Adonis_T1T2_",cohort,".txt")
  adonis.file.name.pdf<-paste0("output/",taxaLevel,"/",taxaLevel,"_Adonis_T1T2_",cohort,".pdf")
  adonis.result<-perform.adonis.all.vars(otu.AN,map.AN,variables,file.Name = adonis.file.name)
  adonis.plot<-plot.adonis(adonis.result,variables.names,adonis.file.name.pdf,show.legend = FALSE)

  #For change in variables
  adonis.file.name<-paste0("output/",taxaLevel,"/",taxaLevel,"_Adonis_T1_",cohort,".txt")
  adonis.file.name.pdf<-paste0("output/",taxaLevel,"/",taxaLevel,"_Adonis_T1_",cohort,".pdf")
  adonis.result<-perform.adonis.all.vars(otu.AN.baseline,map.AN.baseline,variables.change,file.Name = adonis.file.name)
  adonis.plot<-plot.adonis(adonis.result,variables.change.names,adonis.file.name.pdf,show.legend = FALSE)
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


