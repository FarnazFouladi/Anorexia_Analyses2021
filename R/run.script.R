rm(list=ls())

set.seed(120)

source("R/load.packages.R")
source("R/load.data.R")
source("R/sample.selection.R")
source("R/plot.R")
source("R/statistics.R")
source("R/save.plot.R")
source("R/run.analysis.R")

taxa<-c("phylum","class","order","family","genus","species")

#path to metadata and count tables
meta.path="data/meta_master.txt"
otu.path="data/anorexia2020Sep17_taxaCount_%s.tsv"
pathway.path="data/humanN2_pathabundance_cpm.tsv"


#******************Analyses on all taxonomic levels and pathways****************
lapply(taxa,function(t) run.analysis.all(t,otu.path,meta.path,is.taxonomy = TRUE))
run.analysis.all("pathway",pathway.path,meta.path,is.taxonomy = FALSE)
