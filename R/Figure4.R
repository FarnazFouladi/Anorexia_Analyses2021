# Figure 4 ----------------------------------------------------------------

dir.name <- 'Figures/Figure4'
if(!dir.exists(dir.name)){
  dir.create(dir.name,recursive = TRUE,showWarnings = FALSE)
}

variables<-c("BMIchangePerDay")
variable.names<-c("BMI change per day (kg/m2)")
t="Species"

getMLMResults<-function(cohort,t){
  df<-lapply(variables,function(x) read.table(paste0("output/",t,"/",t,"_MLM_",x,"_",cohort,".txt"),sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote = ""))

  names(df)<-variables
  df1<-lapply(df,function(x){
    x$adjusted.p.variable.transformed = getlog10p(x$`adjusted.p-variable`,x$slop)
    return(x)})


  df.combined<-data.frame(taxa=rownames(df1[[1]]))
  for (i in 1:length(df)){
    df.combined<-cbind(df.combined,df1[[i]]$adjusted.p.variable.transformed)
  }
  df.combined<-df.combined %>% tibble::column_to_rownames("taxa")
  colnames(df.combined)<-paste0(variable.names,"_",ifelse(cohort == 'All','All',ifelse(cohort == 'UNC','CEED','ACUTE')))
  return(df.combined)
}

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


pdf(paste0("Figures/Figure4/",t,"_Heatmap_MLM_most_abun.pdf"),height = 15,width = 10)
Heatmap(sig.result1,bottom_annotation = ha,row_names_gp = gpar(fontsize = 10),
        col=col_fun,heatmap_height = unit(25, "cm"),heatmap_width = unit(12, "cm"),
        name = "log10 p-value",cluster_columns = FALSE)

dev.off()

#Scatterplots
t="Species"
cohort="Denver"

df<-read.table(paste0("output/",t,"/",t,"_MLM_BMIchangePerDay_",cohort,".txt"),sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote = "")
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

pdf(paste0("Figures/Figure4/",t,"_ScatterPlots_MLM.pdf"),height = 10,width = 7)
grid.arrange(plots.genus[[1]],plots.genus[[2]],plots.genus[[3]],
             plots.genus[[4]],plots.genus[[5]],plots.genus[[6]],nrow=3,ncol=2)
dev.off()
