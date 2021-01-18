#Figures for the manuscript


# Prepare data ------------------------------------------------------------

meta.path="data/meta_master.txt"
taxa.table.genus="data/anorexia2020Sep17_taxaCount_norm_Log10_genus.tsv"
taxa.table.species="data/anorexia2020Sep17_taxaCount_norm_Log10_species.tsv"
pathway.table="data/humanN2_pathabundance_cpm.tsv"


genus<-load.data(meta.path,taxa.table.genus,normalize = FALSE)[[2]]
genus.unc<-select.samples(map,genus,site = "UNC")[[2]]
genus.denver<-select.samples(map,genus,site = "Denver")[[2]]
species<-load.data(meta.path,taxa.table.species,normalize = FALSE)[[2]]
species.unc<-select.samples(map,species,site = "UNC")[[2]]
species.denver<-select.samples(map,species,site = "Denver")[[2]]

map=load.data(meta.path,taxa.table.genus,normalize = FALSE)[[1]]
map.unc<-select.samples(map,genus,site = "UNC")[[1]]
map.denver<-select.samples(map,genus,site = "Denver")[[1]]

pathway<-load.pathways(meta.path,pathway.table)[[2]]
pathway.unc<-select.samples(map,species,site = "UNC")[[2]]
pathway.denver<-select.samples(map,species,site = "Denver")[[2]]


# Figure1 -----------------------------------------------------------------

plot1<-pco.plot(map,species,color = "Type", shape="Location" ,plot.title = "Species",show.legend = FALSE)
plot2<-pco.plot(map.denver,species.denver,color = "Type", shape="Location" ,plot.title = "Species, Denver",show.legend = FALSE)
plot3<-pco.plot(map.unc,species.unc,color = "Type", shape="Location" ,plot.title = "Species, UNC",show.legend = FALSE)
plot4<-pco.plot(map,pathway,color = "Type", shape="Location" ,plot.title = "Pathways",show.legend =FALSE)
plot5<-pco.plot(map.denver,pathway.denver,color = "Type", shape="Location" ,plot.title = "Pathways, Denver",show.legend = FALSE)
plot6<-pco.plot(map.unc,pathway.unc,color = "Type", shape="Location" ,plot.title = "Pathways, UNC",show.legend = FALSE)

pdf("./output/manuscript_figures/Figure1.pdf",width = 10,height = 7)
grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6,nrow=2,ncol=3)
dev.off()


# Figure 2 ----------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)

t="Phylum"

t.test.result<-read.table(paste0("output/",t,"_t-test_All.txt"),sep = "\t",header = TRUE,row.names = 1,check.names = FALSE,quote = "")
t.test.result.unc<-read.table(paste0("output/",t,"_t-test_UNC.txt"),sep = "\t",header = TRUE,row.names = 1,check.names = FALSE,quote = "")
t.test.result.denver<-read.table(paste0("output/",t,"_t-test_Denver.txt"),sep = "\t",header = TRUE,row.names = 1,check.names = FALSE,quote = "")

result.all<-cbind(t.test.result,t.test.result.unc,t.test.result.denver)
colnames(result.all)<-c(paste0(colnames(t.test.result),".All"),paste0(colnames(t.test.result.unc),".UNC"),
                        paste0(colnames(t.test.result.denver),".Denver"))

pvals<-result.all %>% select(adj.p.vals.HC.T1.All,adj.p.vals.HC.T2.All,adj.p.vals.T1.T2.All,
                                  adj.p.vals.HC.T1.UNC,adj.p.vals.HC.T2.UNC,adj.p.vals.T1.T2.UNC,
                                  adj.p.vals.HC.T1.Denver,adj.p.vals.HC.T2.Denver,adj.p.vals.T1.T2.Denver)

t.stats<-result.all %>% select(t.HC.T1.All,t.HC.T2.All,t.T1.T2.All,
                               t.HC.T1.UNC,t.HC.T2.UNC,t.T1.T2.UNC,
                               t.HC.T1.Denver,t.HC.T2.Denver,t.T1.T2.Denver)


logpvals<-sapply(1:9,function(x){pvals[,x]<- -log10(pvals[,x])})
colnames(logpvals)<-c("HC vs. T1 - All","HC vs. T2 - All","T1 vs. T2 - All",
                      "HC vs. T1 - UNC","HC vs. T2 - UNC","T1 vs. T2 - UNC",
                      "HC vs. T1 - Denver","HC vs. T2 - Denver","T1 vs. T2 - Denver")
rownames(logpvals)<-rownames(pvals)

#sig.result<-logpvals[rowSums(logpvals > 1.30103 )>0,]

cohorts<-factor(c(rep("All",3),rep("UNC",3),rep("Denver",3)))
cohort.cols = list(Cohorts = c("All" = "yellow", "UNC" = "purple","Denver" = "green"))

ha <- HeatmapAnnotation(
  Cohort = cohorts, col = cohort.cols
)

p.cols=brewer.pal(5,"OrRd")
col_fun = colorRamp2(c(1.3,2,3,4,5), p.cols)

pdf(paste0("output/manuscript_figures/",t,"_Heatmap.pdf"),height = 12,width = 10)
Heatmap(logpvals,col=col_fun,top_annotation = ha,row_names_gp = gpar(fontsize = 10),heatmap_height = unit(30, "cm"),
        name = "log10 p-value")

dev.off()


