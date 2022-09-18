# Figure 3 ----------------------------------------------------------------
#Compare p-values between cohorts:

source("../resources/PrepareDataForFigures.R")

moduleDir = dirname(getwd())
pipeRoot = dirname(moduleDir)
outDir = file.path(moduleDir, "output")
message("Output directory: ", outDir)

# Upstream modules
# taxaModuleName = dir(pipeRoot, pattern = "TaxonomyAnalysis", full.names = FALSE)
# pathwayModuleName = dir(pipeRoot, pattern = "PathwayAnalysis", full.names = FALSE)

if(!dir.exists(outDir)){
  dir.create(outDir,recursive = TRUE,showWarnings = FALSE)
}

plot.list <- list()
index <- 1

for(t in c("Genus","Pathway")){

  site<- c("UNC","Denver")
  site.names <- c("CEED","ACUTE")
  df.groups<-lapply(site, function(x) load.ttest.results(t, paste0("_t-test_",x,".txt")) )
  df.groups1<-lapply(df.groups, function(x){
    x$p.vals.HC.T1.transformed <- getlog10p(x$p.vals.HC.T1,x$t.HC.T1)
    x$p.vals.HC.T2.transformed <- getlog10p(x$p.vals.HC.T2,x$t.HC.T2)
    x$p.vals.T1.T2.transformed <- getlog10p(x$p.vals.T1.T2,x$t.T1.T2)
    return(x)
  })

  plot1<-plot.pvals(df.groups1[[1]]$p.vals.HC.T1.transformed,
                    df.groups1[[2]]$p.vals.HC.T1.transformed,
                    site.names[1],site.names[2],plot.title = "non-ED vs T1")

  plot2<-plot.pvals(df.groups1[[1]]$p.vals.HC.T2.transformed,
                    df.groups1[[2]]$p.vals.HC.T2.transformed,
                    site.names[1],site.names[2],plot.title = "non-ED vs T2")

  plot3<-plot.pvals(df.groups1[[1]]$p.vals.T1.T2.transformed,
                    df.groups1[[2]]$p.vals.T1.T2.transformed,
                    site.names[1],site.names[2],plot.title = "T1 vs T2")

  plot.list[[index]] <- plot1
  index <- index + 1
  plot.list[[index]] <- plot2
  index <- index + 1
  plot.list[[index]] <- plot3
  index <- index + 1


}

outFile=file.path(outDir, "Figure3.pdf")
pdf(outFile,width = 12,height = 8)
plot_grid(plotlist = plot.list,nrow = 2, ncol = 3)
dev.off()
