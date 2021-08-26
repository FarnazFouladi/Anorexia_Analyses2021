#Figures for the manuscript

set.seed(120)

source("../resources/load.packages.R")
source("../resources/load.data.R")
source("../resources/sample.selection.R")
source("../resources/statistics.R")
source("../resources/plot.R")
source("../resources/colors.R")

# Locate data ------------------------------------------------------------
moduleDir = dirname(getwd())
pipeRoot = dirname(moduleDir)
meta.path = file.path(pipeRoot, "meta_master.txt")

taxa.table.order = file.path(pipeRoot, "input","taxaCounts", "anorexia2020Sep17_taxaCount_order.tsv")
taxa.table.genus = file.path(pipeRoot, "input","taxaCounts", "anorexia2020Sep17_taxaCount_genus.tsv")
taxa.table.species = file.path(pipeRoot, "input","taxaCounts", "anorexia2020Sep17_taxaCount_species.tsv")
pathway.table = file.path(pipeRoot, "input","pathway", "humanN2_pathabundance_cpm.tsv")


# Prepare data ------------------------------------------------------------
order<-load.data(meta.path,taxa.table.order,normalize = TRUE)[[2]]
genus<-load.data(meta.path,taxa.table.genus,normalize = TRUE)[[2]]
map=load.data(meta.path,taxa.table.genus,normalize = TRUE)[[1]]
genus.unc<-select.samples(map,genus,site = "UNC")[[2]]
genus.denver<-select.samples(map,genus,site = "Denver")[[2]]
species<-load.data(meta.path,taxa.table.species,normalize = TRUE)[[2]]
species.unc<-select.samples(map,species,site = "UNC")[[2]]
species.denver<-select.samples(map,species,site = "Denver")[[2]]

#load unnormalized tables to calculate shannon diversity
species_no_normalized<-load.data(meta.path,taxa.table.species,normalize = FALSE)[[2]]
map$shannon_taxonomy <- diversity(species_no_normalized)

pathway<-load.pathways(meta.path,pathway.table)[[2]]
pathway.unc<-select.samples(map,pathway,site = "UNC")[[2]]
pathway.denver<-select.samples(map,pathway,site = "Denver")[[2]]

map$shannon_pathways <- diversity(pathway)
map.unc<-select.samples(map,genus,site = "UNC")[[1]]
map.denver<-select.samples(map,genus,site = "Denver")[[1]]

##Subsetting the tables

#Anorexia patients
map.AN<-map %>% filter(Type != "HC")
species.AN<-species %>% filter(map$Type != "HC")
genus.AN<-genus %>% filter(map$Type != "HC")
pathway.AN<-pathway %>% filter(map$Type != "HC")
pathway.AN.baseline <- pathway.AN %>% filter(map.AN$Type == "T1")


#Anorexia baseline
map.AN.baseline <- map.AN %>% filter(Type == "T1")
species.AN.baseline <- species.AN %>% filter(map.AN$Type == "T1")
genus.AN.baseline <- genus.AN %>% filter(map.AN$Type == "T1")

#Anorexia baseline at UNC and Denver
map.AN.baseline.UNC <-map.AN.baseline %>% filter(Location=="UNC")
map.AN.baseline.Denver <-map.AN.baseline %>% filter(Location=="Denver")

#variables
unc.meta.variables<-c("BMI","TotalFatPercent","StandardBIA_FatPercent",
                      "Trunk..Percent.Fat","Head..Percent.Fat",
                      "Trunk..Fat.mass..g.","Head..Fat.mass..g.","Total..Fat.mass..g.")

unc.meta.variables.names<-c("BMI (kg/m2)","Total Fat (percent)","Standard BIA Fat (percent)",
                            "Trunk Fat (percent)","Head Fat (percent)",
                            "Trunk Fat Mass (g)","Head Fat Mass (g)","Total Fat Mass (g)")

unc.meta.change<-c("BMIchange","BMIchangePerDay","trunkFatChange","headFatChange",
                   "totalFatChange","trunkPercentFatChange","headPercentFatChange",
                   "totalPercentFatChange")

unc.meta.change.names<-c("BMI Change (kg/m2)","BMI Change per Day (kg/m2.day)",
                         "Trunk Fat Change (g)","Head Fat Change (g)", "Total Fat Change (g)",
                         "Trunk Fat Change (percent)","Head Fat Change (percent)", "Total Fat Change (percent)")

denver.meta.change<-c("BMIchange","BMIchangePerDay")
denver.meta.change.names<-c("BMI Change (kg/m2)","BMI Change per Day (kg/m2.day)")

