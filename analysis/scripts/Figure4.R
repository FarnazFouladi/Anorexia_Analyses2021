# Figure 4 ----------------------------------------------------------------

source("../resources/PrepareDataForFigures.R")

moduleDir = dirname(getwd())
pipeRoot = dirname(moduleDir)
outDir = file.path(moduleDir, "output")
message("Output directory: ", outDir)

# Upstream modules
taxaModuleName = dir(pipeRoot, pattern = "TaxonomyAnalysis", full.names = FALSE)
# pathwayModuleName = dir(pipeRoot, pattern = "PathwayAnalysis", full.names = FALSE)

if(!dir.exists(outDir)){
  dir.create(outDir,recursive = TRUE,showWarnings = FALSE)
}

variables<-c("BMIchangePerDay")
variable.names<-c("BMI change per day (kg/m2)")
t="Species"

getMLMResults<-function(cohort,t,transform.pvals = TRUE){
  df<-lapply(variables,function(x) load.MLM.results(t, paste0("_MLM_",x,"_",cohort,".txt")))

  names(df)<-variables

  if(transform.pvals)
    df<-lapply(df,function(x){
      x$adjusted.p.variable.transformed = getlog10p(x$`adjusted.p-variable`,x$slop)
      return(x)})

  df.combined<-data.frame(taxa=rownames(df[[1]]))
  for (i in 1:length(df)){
    df.combined<-cbind(df.combined,df[[i]][,ncol(df[[i]])])
  }
  df.combined<-df.combined %>% tibble::column_to_rownames("taxa")
  colnames(df.combined)<-paste0(variable.names,"_",ifelse(cohort == 'All','All',ifelse(cohort == 'UNC','CEED','ACUTE')))
  return(df.combined)
}

#Get pvalues and write the table
df.all<-getMLMResults("All",t,transform.pvals = FALSE)
df.unc<-getMLMResults("UNC",t,transform.pvals = FALSE)
df.denver<-getMLMResults("Denver",t,transform.pvals = FALSE)

df.final_non_transformed<-cbind(df.all,df.unc,df.denver)
df.final_non_transformed <- df.final_non_transformed %>% tibble::rownames_to_column("Taxa")
MLM.table.file=file.path(outDir, paste0(t,"_MLM_resulst.txt"))
write.table(df.final_non_transformed, MLM.table.file, sep = '\t',quote=FALSE,row.names = FALSE)

#Get pvalues, log10 transform, and plot a heatmap
df.all<-getMLMResults("All",t)
df.unc<-getMLMResults("UNC",t)
df.denver<-getMLMResults("Denver",t)

df.final<-cbind(df.all,df.unc,df.denver)

sig.result<-as.matrix(df.final[rowSums(df.final > -log10(0.05) | df.final < log10(0.05) )>0,])

#Keep the first 50 most abundant taxa for Denver:
cohort = 'Denver'
otu.AN.baseline <- get(paste0(tolower(t),'.AN.baseline'))

map.AN.baseline.cohort<-map.AN.baseline %>% filter(Location==cohort)
otu.AN.baseline.cohort<-otu.AN.baseline %>% filter(map.AN.baseline$Location==cohort)
otu_sig<-otu.AN.baseline.cohort[,rownames(sig.result)]
otu_ordered<-otu_sig[,order(colSums(otu_sig),decreasing = TRUE)]
sig.result1<-sig.result[match(colnames(otu_ordered),rownames(sig.result)),]
sig.result1<-sig.result1[1:50,]


cohorts<-factor(c(rep("All",length(variables)),rep("CEED",length(variables)),rep("ACUTE",length(variables))),levels = c('All','CEED','ACUTE'))
cohort.cols = list(Cohorts = c("All" = "#FFD92F", "CEED" = "#F0027F","ACUTE" = "#1F78B4"))

ha <- HeatmapAnnotation(
  Cohorts = cohorts, col = cohort.cols
)

p.cols.red=brewer.pal(3,"OrRd")
pcols.blue=brewer.pal(3,"Blues")
col_fun = colorRamp2(c(-1.5,-1,-0.5,0.5,1,1.5),c("blue","#9ECAE1","#DEEBF7","#FEE8C8","#FDBB84","red"))

heatmap.file=file.path(outDir, paste0(t,"_Heatmap_MLM_most_abun.pdf"))
pdf(heatmap.file,height = 15,width = 10)
Heatmap(sig.result1,bottom_annotation = ha,row_names_gp = gpar(fontsize = 10),
        col=col_fun,heatmap_height = unit(25, "cm"),heatmap_width = unit(12, "cm"),
        name = "log10 p-value",cluster_columns = FALSE)

dev.off()

#Scatterplots
t="Species"
cohort="Denver"

mlm.results=file.path(pipeRoot, taxaModuleName, "output", t, paste0(t,"_MLM_BMIchangePerDay_",cohort,".txt"))
message("Reading file: ", mlm.results)
df<-read.table(mlm.results, sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote = "")
df<-df[df$`adjusted.p-variable`<0.05,]
otu_sig<-otu.AN.baseline.cohort[,rownames(df)]
otu_ordered<-otu_sig[,order(colSums(otu_sig),decreasing = TRUE)]
df1<-df[match(colnames(otu_ordered),rownames(df)),]


taxa <- c('Intestinimonas butyriciproducens','Bacteroides sp. A1C1','Pseudomonas stutzeri',
          'Clostridium perfringens','[Ruminococcus] gnavus','Lachnospiraceae bacterium')

df2 <- df1[taxa,]
pvals <- df2$`adjusted.p-variable`
xlab="BMI change per Day (kg/m2)"
plots.genus<-lapply(1:6, function(x) scatter.plot(map.AN.baseline.cohort,otu.AN.baseline.cohort,taxa = taxa[x],"BMIchangePerDay",xlab,pvals[x]))

scatter.file=file.path(outDir, paste0(t,"_ScatterPlots_MLM.pdf"))
pdf(scatter.file,height = 10,width = 7)
grid.arrange(plots.genus[[1]],plots.genus[[2]],plots.genus[[3]],
             plots.genus[[4]],plots.genus[[5]],plots.genus[[6]],nrow=3,ncol=2)
dev.off()
