# Figure 2 ----------------------------------------------------------------

dir.name <- 'Figures/Figure2'
if(!dir.exists(dir.name)){
  dir.create(dir.name,recursive = TRUE,showWarnings = FALSE)
}

##########Heatmaps##########

level = 2

for (t in c("Genus","Species","Pathway")){

  t.test.result<-read.table(paste0("output/",t,"/",t,"_t-test_All.txt"),sep = "\t",header = TRUE,row.names = 1,check.names = FALSE,quote = "")
  t.test.result.unc<-read.table(paste0("output/",t,"/",t,"_t-test_UNC.txt"),sep = "\t",header = TRUE,row.names = 1,check.names = FALSE,quote = "")
  t.test.result.denver<-read.table(paste0("output/",t,"/",t,"_t-test_Denver.txt"),sep = "\t",header = TRUE,row.names = 1,check.names = FALSE,quote = "")

  result.all<-cbind(t.test.result,t.test.result.unc,t.test.result.denver)
  colnames(result.all)<-c(paste0(colnames(t.test.result),".All"),paste0(colnames(t.test.result.unc),".UNC"),
                          paste0(colnames(t.test.result.denver),".Denver"))



  result.all.write <- result.all %>% select(-c("taxa.All" ,"taxa.UNC","taxa.Denver")) %>%
    tibble::rownames_to_column("Taxa")
  write.table(result.all.write, paste0("Figures/Figure2/",t,"_ttest_resulst.txt"), sep = '\t',quote=FALSE,row.names = FALSE)

  pvals<-result.all %>% select(adj.p.vals.HC.T1.All,adj.p.vals.HC.T2.All,adj.p.vals.T1.T2.All,
                               adj.p.vals.HC.T1.UNC,adj.p.vals.HC.T2.UNC,adj.p.vals.T1.T2.UNC,
                               adj.p.vals.HC.T1.Denver,adj.p.vals.HC.T2.Denver,adj.p.vals.T1.T2.Denver)

  t.stats<-result.all %>% select(t.HC.T1.All,t.HC.T2.All,t.T1.T2.All,
                                 t.HC.T1.UNC,t.HC.T2.UNC,t.T1.T2.UNC,
                                 t.HC.T1.Denver,t.HC.T2.Denver,t.T1.T2.Denver)


  logpvals<-sapply(1:9,function(x){pvals[,x]<- -log10(pvals[,x])})
  colnames(logpvals)<-c("non-ED vs. T1 - All","non-ED vs. T2 - All","T1 vs. T2 - All",
                        "non-ED vs. T1 - CEED","non-ED vs. T2 - CEED","T1 vs. T2 - CEDD",
                        "non-ED vs. T1 - ACUTE","non-ED vs. T2 - ACUTE","T1 vs. T2 - ACUTE")
  rownames(logpvals)<-rownames(pvals)

  sig.result<-logpvals[rowSums(logpvals > -log10(0.05) )>0,]


  #Keep the first 50 most abundant taxa:
  otu <- get(tolower(t))
  otu<-otu[,rownames(sig.result)]
  ordered_otu<-otu[,order(colMeans(otu),decreasing = TRUE)]
  sig.result1<-sig.result[match(colnames(ordered_otu),rownames(sig.result)),]
  sig.result1<-sig.result1[1:50,]

  cohorts<-factor(c(rep("All",3),rep("CEED",3),rep("ACUTE",3)),levels = c('All','CEED','ACUTE'))
  cohort.cols = list(Cohorts = c("All" = "#FFD92F", "CEED" = "#F0027F","ACUTE" = "#1F78B4"))

  ha <- HeatmapAnnotation(
    Cohorts = cohorts, col = cohort.cols
  )

  p.cols=brewer.pal(5,"OrRd")
  col_fun = colorRamp2(c(1.3,2,3,4,5), p.cols)


  #Add annontaion for pathways
  if(t == 'Pathway'){
    annotated_path<- readxl::read_xlsx('output/Pathway/significant_pathways.xlsx',na = 'NA')
    annotated_path$level1 <- sapply(as.character(annotated_path$Class),function(x) {if (x!='Superpathways') strsplit(x,';')[[1]][1] else return('Superpathways')})
    annotated_path$level2 <- sapply(as.character(annotated_path$Class),function(x) {if (x!='Superpathways') strsplit(x,';')[[1]][2] else return('Superpathways')})

    #Select pathways that we want ti show in the heat map (the top 50 abundant pathways that were significant)
    annotated_path_heatmap <- annotated_path %>% filter(Pathways %in% rownames(sig.result1))
    annotated_path_heatmap <- annotated_path_heatmap[match(rownames(sig.result1),annotated_path_heatmap$Pathways),]

    rownames(sig.result1) <- sapply(rownames(sig.result1),function(x){strsplit(x,': ')[[1]][2]})

    if(level == 1){
      pdf(paste0("Figures/Figure2/",t,"_Heatmap_level1.pdf"),height = 12,width = 17)

      row_ha = rowAnnotation(Pathways = annotated_path_heatmap$level1,
                             col=list(Pathways = c("Biosynthesis" = "#66C2A5",
                                                   "Generation of Precursor Metabolites and Energy"="#FC8D62",
                                                   "Degradation/Utilization/Assimilation" = "#BC80BD",
                                                   "Superpathways" = "#5d8fd5" )),annotation_name_side = "top")
    } else{
      pdf(paste0("Figures/Figure2/",t,"_Heatmap_level2.pdf"),height = 12,width = 17)
      #level2
      myCols <- c("#b4426b","#d395a5","#de424e","#973a35","#855d55","#e08f71",
                  "#ded43e","#63684b","#cde18a","#b4e445","#74a737","#c8d6bd",
                  "#66db4a","#6cdbaf","#5b9c86","#3b6761","#bec3dc","#59559b",
                  "#e4408f","#d04ad1","#943d93","#b89ad8","#783fcf","#6dd5db",
                  "#3a6ce3","#5d8fd5")

      level2.cols <- myCols[1:length(unique(annotated_path_heatmap$level2))]
      names(level2.cols) <- unique(annotated_path_heatmap$level2)

      row_ha = rowAnnotation(Pathways = annotated_path_heatmap$level2, col=list(Pathways=level2.cols),
                             annotation_name_side = "top")
    }


    hm<-Heatmap(sig.result1,col=col_fun,bottom_annotation = ha,
                right_annotation =  row_ha,row_names_gp = gpar(fontsize = 10),
                heatmap_height = unit(25, "cm"),heatmap_width = unit(15, "cm"),
                name = "log10 p-value",cluster_columns = FALSE)

    draw(hm,heatmap_legend_side = "left",annotation_legend_side = "left")
  } else{

    pdf(paste0("Figures/Figure2/",t,"_Heatmap.pdf"),height = 12,width = 17)

    hm<-Heatmap(sig.result1,col=col_fun,bottom_annotation = ha,
                row_names_gp = gpar(fontsize = 10),heatmap_height = unit(25, "cm"),
                heatmap_width = unit(15, "cm"),
                name = "log10 p-value",cluster_columns = FALSE)

    draw(hm,heatmap_legend_side = "left",annotation_legend_side = "left")


  }
  dev.off()
}

# Extended Figure 2 and 3 ------------------------------------------------------
#Boxplots for the most significant species and pathways in the heatmap

taxa.df <- cbind(species,map)
taxa.df.denver <- cbind(species.denver,map.denver)
taxa.df.unc <- cbind(species.unc,map.unc)

pathway.df <- cbind(pathway,map)
pathway.df.denver <- cbind(pathway.denver,map.denver)
pathway.df.unc <- cbind(pathway.unc,map.unc)


plots.all <- list()
index <- 1

for (t in c("Species","Pathway")){

  t.test.result<-read.table(paste0("output/",t,"/",t,"_t-test_All.txt"),sep = "\t",header = TRUE,row.names = 1,check.names = FALSE,quote = "")
  t.test.result.unc<-read.table(paste0("output/",t,"/",t,"_t-test_UNC.txt"),sep = "\t",header = TRUE,row.names = 1,check.names = FALSE,quote = "")
  t.test.result.denver<-read.table(paste0("output/",t,"/",t,"_t-test_Denver.txt"),sep = "\t",header = TRUE,row.names = 1,check.names = FALSE,quote = "")

  if (t == "Species"){
    names <- c("Bifidobacterium adolescentis","Flavonifractor plautii","Collinsella aerofaciens","Faecalibacterium prausnitzii")
    is.taxonomy = TRUE
    data.all <- taxa.df
    data.all.denver <- taxa.df.denver
    data.all.unc <- taxa.df.unc

  }else{
    names<-c("dTDP-L-rhamnose biosynthesis I",
             "Calvin-Benson-Bassham cycle",
             "sucrose degradation III (sucrose invertase)",
             "tRNA charging")
    is.taxonomy = FALSE
    data.all <- pathway.df
    data.all.denver <- pathway.df.denver
    data.all.unc <- pathway.df.unc

    #Remove unmapped reads
    t.test.result <- t.test.result[-c(1,2),]
    t.test.result.denver <- t.test.result[-c(1,2),]
    t.test.result.unc <- t.test.result[-c(1,2),]

    #Trim pathway names
    t.test.result$taxa <- sapply(t.test.result$taxa,function(x){strsplit(x,': ')[[1]][2]})
    t.test.result.denver$taxa <- sapply(t.test.result.denver$taxa,function(x){strsplit(x,': ')[[1]][2]})
    t.test.result.unc$taxa <- sapply(t.test.result.unc$taxa,function(x){strsplit(x,': ')[[1]][2]})

    rownames(t.test.result)<-t.test.result$taxa
    rownames(t.test.result.denver)<-t.test.result.denver$taxa
    rownames(t.test.result.unc)<-t.test.result.unc$taxa

    #Trim pathway names
    colnames(data.all)[1:326]<-sapply(colnames(data.all)[1:326],function(x){strsplit(x,': ')[[1]][2]})
    colnames(data.all.denver)[1:326]<-sapply(colnames(data.all.denver)[1:326],function(x){strsplit(x,': ')[[1]][2]})
    colnames(data.all.unc)[1:326]<-sapply(colnames(data.all.unc)[1:326],function(x){strsplit(x,': ')[[1]][2]})
  }


  #UNC and Denver
  t.test.result.plot <- t.test.result[names,]
  plots<-apply(t.test.result.plot,1,function(x) box.plot(data.all,x,show.legend=FALSE,is.taxonomy))
  for (i in 1:4){
    plots.all[[index]] <- plots[[i]]
    index <- index + 1
  }

  #Denver
  t.test.result.plot.denver <- t.test.result.denver[names,]
  plots.denver<-apply(t.test.result.plot.denver,1,function(x) box.plot(data.all.denver,x,show.legend=FALSE,is.taxonomy))
  for (i in 1:4){
    plots.all[[index]] <- plots.denver[[i]]
    index <- index + 1
  }

  #UNC
  t.test.result.plot.unc <- t.test.result.unc[names,]
  plots.unc<-apply(t.test.result.plot.unc,1,function(x) box.plot(data.all.unc,x,show.legend=FALSE,is.taxonomy))
  for (i in 1:4){
    plots.all[[index]] <- plots.unc[[i]]
    index <- index + 1
  }
}


pdf("Figures/Figure2/Extended_Figure2.pdf",height = 13,width = 14)
gridExtra::grid.arrange(plots.all[[1]],plots.all[[5]],plots.all[[9]],
                   plots.all[[2]],plots.all[[6]],plots.all[[10]],
                   plots.all[[3]],plots.all[[7]],plots.all[[11]],
                   plots.all[[4]],plots.all[[8]],plots.all[[12]],ncol=3)
dev.off()

pdf("Figures/Figure2/Extended_Figure3.pdf",height = 13,width = 14)
gridExtra::grid.arrange(plots.all[[13]],plots.all[[17]],plots.all[[21]],
                        plots.all[[14]],plots.all[[18]],plots.all[[22]],
                        plots.all[[15]],plots.all[[19]],plots.all[[23]],
                        plots.all[[16]],plots.all[[20]],plots.all[[24]],ncol=3)
dev.off()








