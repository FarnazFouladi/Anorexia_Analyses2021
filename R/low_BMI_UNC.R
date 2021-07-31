#Correlation between change in BMI and microbiome for UNC patients with low BMI
rm(list=ls())

source("R/load.data.R")
source("R/sample.selection.R")
source("R/plot.R")
source("R/statistics.R")
source("R/save.plot.R")
source("R/run.analysis.R")

taxa<-c("phylum","class","order","family","genus","species")
t="species"
is.taxonomy = TRUE
cohort = 'UNC'

#path to metadata and count tables
meta.path="data/meta_master.txt"
otu.path="data/anorexia2020Sep17_taxaCount_%s.tsv"
pathway.path="data/humanN2_pathabundance_cpm.tsv"

taxaLevel<-capitalize(t)

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


map.unc<-select.samples(map,otu,site = "UNC")[[1]]
otu.unc<-select.samples(map,otu,site = "UNC")[[2]]
map<-map.unc
otu<-otu.unc

variables.change<-c("BMIchange","BMIchangePerDay",
                    "trunkFatChange","headFatChange",
                    "totalFatChange","trunkPercentFatChange","headPercentFatChange",
                    "totalPercentFatChange")
variables.change.names<-c("BMI Change (kg/m2)","BMI change per Day (kg/m2)",
                          "Trunk Fat Change (g)","Head Fat Change (g)", "Total Fat Change (g)",
                          "Trunk Fat Change (percent)","Head Fat Change (percent)", "Total Fat Change (percent)")

title=paste0(taxaLevel,"_",cohort)

# Baseline Microbiome -----------------------------------------------------
map.AN<-map %>% filter(Type != "HC")
otu.AN<-otu %>% filter(map$Type != "HC")

map.AN.baseline <- map.AN %>% filter(Type == "T1")
otu.AN.baseline <- otu.AN %>% filter(map.AN$Type == "T1")

#Select low BMI
map.AN.baseline.low <- map.AN.baseline %>% filter(BMI < 15)
otu.AN.baseline.low <- otu.AN.baseline %>% filter(map.AN.baseline$BMI < 15)
otu.AN.baseline.low <- otu.AN.baseline.low[,colSums(otu.AN.baseline.low)!=0]


# Mixed Linear Regression Analyses ----------------------------------------

## For change in variables including only T1 microbiome
MLM.file.name<-paste0("output/",taxaLevel,"/",taxaLevel,"_MLM_",variables.change,"_",cohort,"_low_BMI.txt")
MLM.file.name.pdf<-paste0("output/",taxaLevel,"/",taxaLevel,"_MLM_",variables.change,"_",cohort,"_low_BMI.pdf")
MLM.result<-perform.MLM.all.vars(map.AN.baseline,otu.AN.baseline,variables.change,MLM.file.name,changeInVariable = TRUE)
MLM.plots<-lapply(1:length(variables.change),function(x) get.scatter.plots(map.AN.baseline,otu.AN.baseline,variables.change[x],
                                                                           variables.change.names[x],result.test = MLM.result[[x]],legend.show = FALSE))

bug = otu.AN.baseline.low[,1]
variable = map.AN.baseline.low[,"BMIchangePerDay"]
cohort = map.AN.baseline.low[,"Location"]
myData <- data.frame(bug,variable ,cohort)
fit <- summary(lme(bug ~ variable, random = ~ 1 | cohort, data=myData,na.action = na.omit))
fit <- summary(lm(bug ~ variable, data=myData,na.action = na.omit))



