source("R/PrepareDataForFigures.R")
# Figure 3 ----------------------------------------------------------------

variables<-c("BMI","BMIChange","BMIchangePerDay","TimeBetweenRecovery")
variable.names<-c("BMI (kg/m2)","BMI change (kg/m2)","BMI change per day (kg/m2)","Recovery (d)")
t="Pathway"

getMLMResults<-function(cohort,t){
  df<-lapply(variables,function(x) read.table(paste0("output/",t,"_MLM_",x,"_",cohort,".txt"),sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote = ""))

  names(df)<-variables
  df1<-lapply(df,function(x){
    x$adjusted.p.variable.transformed = getlog10p(x$`adjusted.p-variable`,x$slop)
    return(x)})

  df.combined<-data.frame(df1[[1]]$adjusted.p.variable.transformed,
                          df1[[2]]$adjusted.p.variable.transformed,
                          df1[[3]]$adjusted.p.variable.transformed,
                          df1[[4]]$adjusted.p.variable.transformed)
  colnames(df.combined)<-paste0(variable.names,"_",cohort)
  rownames(df.combined)<-rownames(df1[[1]])
  return(df.combined)
}

df.all<-getMLMResults("All",t)
df.unc<-getMLMResults("UNC",t)
df.denver<-getMLMResults("Denver",t)

df.final<-cbind(df.all,df.unc,df.denver)

sig.result<-as.matrix(df.final[rowSums(df.final > -log10(0.1) | df.final < log10(0.1) )>0,])
sig.result<-sig.result[rownames(sig.result)!="NA",]

cohorts<-c(rep("All",length(variables)),rep("UNC",length(variables)),rep("Denver",length(variables)))
cohort.cols = list(Cohorts = c("All" = "yellow", "UNC" = "purple","Denver" = "green"))

ha <- HeatmapAnnotation(
  Cohorts = cohorts, col = cohort.cols
)

p.cols.red=brewer.pal(3,"OrRd")
pcols.blue=brewer.pal(3,"Blues")
col_fun = colorRamp2(c(-5,-3,-1.3,1.3,3,5),c(rev(p.cols.red),pcols.blue) )


pdf(paste0("output/manuscript_figures/",t,"_Heatmap_MLM.pdf"),height = 10,width = 10)
Heatmap(sig.result,top_annotation = ha,row_names_gp = gpar(fontsize = 10),heatmap_height = unit(10, "cm"),
        name = "log10 p-value")

dev.off()

#Scatterplots
t="Genus"
df<-read.table(paste0("output/",t,"_MLM_BMIchangePerDay_",cohort,".txt"),sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote = "")
df<-df[order(df$`adjusted.p-variable`),]
pvals<-df$`adjusted.p-variable`[1:6]
taxa<-rownames(df)[1:6]
xlab="BMI change per Day (kg/m2)"
plots.genus<-lapply(1:6, function(x) scatter.plot(map,genus,taxa = taxa[x],"BMIchangePerDay",xlab,pvals[x]))

pdf(paste0("output/manuscript_figures/",t,"_ScatterPlots_MLM.pdf"),height = 10,width = 7)
grid.arrange(plots.genus[[1]],plots.genus[[2]],plots.genus[[3]],
             plots.genus[[4]],plots.genus[[5]],plots.genus[[6]],nrow=3,ncol=2)
dev.off()

t="Pathway"
df<-read.table(paste0("output/",t,"_MLM_BMIchangePerDay_",cohort,".txt"),sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote = "")
df<-df[order(df$`adjusted.p-variable`),]
pvals<-df$`adjusted.p-variable`[1]
taxa<-rownames(df)[1]
xlab="BMI change per Day (kg/m2)"
plots.pathways<-scatter.plot(map,pathway,taxa = taxa,"BMIchangePerDay",xlab,pvals)



