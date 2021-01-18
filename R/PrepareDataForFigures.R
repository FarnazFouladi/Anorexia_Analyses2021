#Figures for the manuscript


# Prepare data ------------------------------------------------------------

meta.path="data/meta_master.txt"
taxa.table.genus="data/anorexia2020Sep17_taxaCount_norm_Log10_genus.tsv"
taxa.table.species="data/anorexia2020Sep17_taxaCount_norm_Log10_species.tsv"
pathway.table="data/humanN2_pathabundance_cpm.tsv"


genus<-load.data(meta.path,taxa.table.genus,normalize = FALSE)[[2]]
map=load.data(meta.path,taxa.table.genus,normalize = FALSE)[[1]]

genus.unc<-select.samples(map,genus,site = "UNC")[[2]]
genus.denver<-select.samples(map,genus,site = "Denver")[[2]]
species<-load.data(meta.path,taxa.table.species,normalize = FALSE)[[2]]
species.unc<-select.samples(map,species,site = "UNC")[[2]]
species.denver<-select.samples(map,species,site = "Denver")[[2]]
map.unc<-select.samples(map,genus,site = "UNC")[[1]]
map.denver<-select.samples(map,genus,site = "Denver")[[1]]

pathway<-load.pathways(meta.path,pathway.table)[[2]]
pathway.unc<-select.samples(map,pathway,site = "UNC")[[2]]
pathway.denver<-select.samples(map,pathway,site = "Denver")[[2]]

