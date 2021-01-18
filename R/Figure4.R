
# Figure 4 ----------------------------------------------------------------

#Compare p-values between cohorts:
t="Genus"

for(t in c("Genus","Pathway")){

  site<-c("UNC","Denver")
  df.groups<-lapply(site, function(x) read.table(paste0("output/",t,"_t-test_",x,".txt"),sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote = ""))
  df.groups1<-lapply(df.groups, function(x){
    x$p.vals.HC.T1.transformed <- getlog10p(x$p.vals.HC.T1,x$t.HC.T1)
    x$p.vals.HC.T2.transformed <- getlog10p(x$p.vals.HC.T2,x$t.HC.T2)
    x$p.vals.T1.T2.transformed <- getlog10p(x$p.vals.T1.T2,x$t.T1.T2)
    return(x)
  })

  plot1<-plot.pvals(df.groups1[[1]]$p.vals.HC.T1.transformed,
                    df.groups1[[2]]$p.vals.HC.T1.transformed,
                    site[1],site[2],plot.title = "HC vs T1")

  plot2<-plot.pvals(df.groups1[[1]]$p.vals.HC.T2.transformed,
                    df.groups1[[2]]$p.vals.HC.T2.transformed,
                    site[1],site[2],plot.title = "HC vs T2")

  plot3<-plot.pvals(df.groups1[[1]]$p.vals.T1.T2.transformed,
                    df.groups1[[2]]$p.vals.T1.T2.transformed,
                    site[1],site[2],plot.title = "T1 vs T2")

  pdf(paste0("output/manuscript_figures/",t,"_pvalueplots_ttest.pdf"),width = 10,height = 4)
  grid.arrange(plot1,plot2,plot3,ncol=3,nrow=1)
  dev.off()
}

#Compare MLM results

  getPlotMLM<-function(t,variable){
    site<-c("UNC","Denver")
    df.groups<-lapply(site, function(x) read.table(paste0("output/",t,"_MLM_",variable,"_",x,".txt"),sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote = ""))
    df.groups1<-lapply(df.groups, function(x){
      x$`p-variable.transformed` <- getlog10p(x$`p-variable`,x$slope)
      return(x)
    })
    plot1<-plot.pvals(df.groups1[[1]]$`p-variable.transformed`,
                      df.groups1[[2]]$`p-variable.transformed`,
                      site[1],site[2],plot.title = variable)
    plot1
  }

  variables<-c("BMI","BMIchange","BMIchangePerDay","TimeBetweenRecovery")
  plots.genus<-lapply(variables,getPlotMLM,t="Genus")
  plots.pathway<-lapply(variables,getPlotMLM,t="Pathway")

  pdf(paste0("output/manuscript_figures/Genus_pvalueplots_MLM.pdf"),width = 10,height = 10)
  grid.arrange(plots.genus[[1]],plots.genus[[2]],
               plots.genus[[3]],plots.genus[[4]])
  dev.off()

  pdf(paste0("output/manuscript_figures/Pathway_pvalueplots_MLM.pdf"),width = 10,height = 10)
  grid.arrange(plots.pathway[[1]],plots.pathway[[2]],
               plots.pathway[[3]],plots.pathway[[4]])
  dev.off()


