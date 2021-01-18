# Figure 2 ----------------------------------------------------------------

library(ComplexHeatmap)
library(circlize)
source("R/PrepareDataForFigures.R")

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

sig.result<-logpvals[rowSums(logpvals > 1.30103 )>0,]

cohorts<-c(rep("All",3),rep("UNC",3),rep("Denver",3))
cohort.cols = list(Cohorts = c("All" = "yellow", "UNC" = "purple","Denver" = "green"))

ha <- HeatmapAnnotation(
  Cohorts = cohorts, col = cohort.cols
)

p.cols=brewer.pal(5,"OrRd")
col_fun = colorRamp2(c(1.3,2,3,4,5), p.cols)

pdf(paste0("output/manuscript_figures/",t,"_Heatmap.pdf"),height = 12,width = 10)
Heatmap(sig.result,col=col_fun,top_annotation = ha,row_names_gp = gpar(fontsize = 10),heatmap_height = unit(30, "cm"),
        name = "log10 p-value")

dev.off()
