# Figure 2 ----------------------------------------------------------------

dir.name <- 'Figures/Figure2'
if(!dir.exists(dir.name)){
  dir.create(dir.name,recursive = TRUE,showWarnings = FALSE)
}

level = 1

for (t in c("Genus","Species","Pathway")){

  t.test.result<-read.table(paste0("output/",t,"/",t,"_t-test_All.txt"),sep = "\t",header = TRUE,row.names = 1,check.names = FALSE,quote = "")
  t.test.result.unc<-read.table(paste0("output/",t,"/",t,"_t-test_UNC.txt"),sep = "\t",header = TRUE,row.names = 1,check.names = FALSE,quote = "")
  t.test.result.denver<-read.table(paste0("output/",t,"/",t,"_t-test_Denver.txt"),sep = "\t",header = TRUE,row.names = 1,check.names = FALSE,quote = "")

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


  #Keep the first 30 most abundant taxa:
  otu <- get(tolower(t))
  otu<-otu[,rownames(sig.result)]
  ordered_otu<-otu[,order(colMeans(otu),decreasing = TRUE)]
  #sig.result1<-sig.result[match(colnames(ordered_otu),rownames(sig.result)),]
  #sig.result1<-sig.result1[1:30,]
  sig.result1 <- sig.result

  cohorts<-c(rep("All",3),rep("UNC",3),rep("Denver",3))
  cohort.cols = list(Cohorts = c("All" = "#FFD92F", "UNC" = "#F0027F","Denver" = "#1F78B4"))

  ha <- HeatmapAnnotation(
    Cohorts = cohorts, col = cohort.cols
  )

  p.cols=brewer.pal(5,"OrRd")
  col_fun = colorRamp2(c(1.3,2,3,4,5), p.cols)


  #Add annontaion for pathways
  if(t == 'Pathway'){
    annotated_path<- readxl::read_xlsx('output/Pathway/significant_pathways.xlsx',na = 'NA')
    annotated_path$level1 <- sapply(as.character(annotated_path$Class),function(x) {strsplit(x,';')[[1]][1]})
    annotated_path$level2 <- sapply(as.character(annotated_path$Class),function(x) {strsplit(x,';')[[1]][2]})
    annotated_path_heatmap <- annotated_path %>% filter(Pathways %in% rownames(sig.result1))
    annotated_path_heatmap <- annotated_path_heatmap[match(rownames(sig.result1),annotated_path_heatmap$Pathways),]

    #Remove NA
    annotated_path_heatmap <- annotated_path_heatmap[!is.na(annotated_path_heatmap$level2),]
    sig.result1 <- sig.result1[annotated_path_heatmap$Pathways,]

    rownames(sig.result1) <- sapply(rownames(sig.result1),function(x){strsplit(x,': ')[[1]][2]})
    #rownames(sig.result1) <- sapply(rownames(sig.result1),function(x){strsplit(x,':')[[1]][2]})
    if(level == 1){
      pdf(paste0("Figures/Figure2/",t,"_Heatmap_level1.pdf"),height = 12,width = 17)
      row_ha = rowAnnotation(Pathways = annotated_path_heatmap$level1, col=list(Pathways = c("Biosynthesis" = "#66C2A5", "Generation of Precursor Metabolites and Energy"="#FC8D62",
                                                                                             "Degradation/Utilization/Assimilation" = "#BC80BD")),annotation_name_side = "top")
    } else{
      pdf(paste0("Figures/Figure2/",t,"_Heatmap_level2.pdf"),height = 12,width = 17)

      #level2

      myCols <- c("#b4426b",
                      "#d395a5",
                      "#de424e",
                      "#973a35",
                      "#855d55",
                      "#e08f71",
                      "#ded43e",
                      "#63684b",
                      "#cde18a",
                      "#b4e445",
                      "#74a737",
                      "#c8d6bd",
                      "#66db4a",
                      "#6cdbaf",
                      "#5b9c86",
                      "#3b6761",
                      "#bec3dc",
                      "#59559b",
                      "#e4408f",
                      "#d04ad1",
                      "#943d93",
                      "#b89ad8",
                      "#783fcf",
                      "#6dd5db",
                      "#3a6ce3",
                      "#5d8fd5",
                      "#FFEA46FF",
                      "#FCFDBFFF")

      pdf(paste0("Figures/Figure2/",t,"_Heatmap_level2_all.pdf"),height = 30,width = 17)


      level2.cols <- myCols[1:length(unique(annotated_path_heatmap$level2))]
      names(level2.cols) <- unique(annotated_path_heatmap$level2)

      row_ha = rowAnnotation(Pathways = annotated_path_heatmap$level2, col=list(Pathways=level2.cols),annotation_name_side = "top")
    }

    hm<-Heatmap(sig.result1,col=col_fun,bottom_annotation = ha,right_annotation =  row_ha,row_names_gp = gpar(fontsize = 10),heatmap_height = unit(65, "cm"),heatmap_width = unit(15, "cm"),
                name = "log10 p-value",cluster_columns = FALSE)

     draw(hm,heatmap_legend_side = "left",annotation_legend_side = "left")
  } else{
    pdf(paste0("Figures/Figure2/",t,"_Heatmap.pdf"),height = 12,width = 17)
    hm<-Heatmap(sig.result1,col=col_fun,bottom_annotation = ha,row_names_gp = gpar(fontsize = 10),heatmap_height = unit(20, "cm"),heatmap_width = unit(15, "cm"),
                name = "log10 p-value",cluster_columns = FALSE)
    draw(hm,heatmap_legend_side = "left",annotation_legend_side = "left")


  }
  dev.off()
}

#
num_pathways <- matrix(c(78,90,14,19),nrow=2,ncol=2,
                       dimnames = list(c('Denver','UNC'),c('HC-T1','HC-T2')))
chisq.test(num_pathways)
fisher_test(num_pathways)
