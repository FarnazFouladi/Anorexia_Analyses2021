rm(list=ls())


source("R/load.data.R")
source("R/sample.selection.R")
source("R/plot.R")
source("R/statistics.R")
source("R/save.plot.R")
source("R/run.analysis.R")



taxa<-c("phylum","order","family","species")

meta.path="data/meta_master.txt"
taxa.table="data/anorexia2020Sep17_taxaCount_norm_Log10_%s.tsv"
pathway.table="data/humanN2_pathabundance_cpm.tsv"


#******************Analyses on all taxonomic levels and pathways****************

taxa="genus"

lapply(taxa,function(x) run.analysis.all(x,taxa.table,meta.path,taxanomy = TRUE))
run.analysis.all("pathway",pathway.table,meta.path,taxanomy = FALSE)





















