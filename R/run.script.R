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

#Set colors and names
cols.type <- brewer.pal(3, "Dark2")
cols.site <- c("#1F78B4","#F0027F")
names.type <- c('non-ED','T1','T2')
names.site <- c("Denver","UNC")

#******************Analyses on all taxonomic levels and pathways****************
lapply(taxa,function(t) run.analysis.all(t,otu.path,meta.path,is.taxonomy = TRUE))
run.analysis.all("pathway",pathway.path,meta.path,is.taxonomy = FALSE)
