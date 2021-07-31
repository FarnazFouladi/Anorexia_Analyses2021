#Ordination plot for figure 1
pco.plot<-function(map,otu,color.group,plot.title,axis1=1,axis2=2,figure.colors,
                   show.legend=TRUE,legend.names=NULL,show.stat=TRUE,show.ellipses=TRUE,
                   pdf.name){

  pco<-capscale(otu ~ 1, distance = "bray")
  percent<-pco$CA$eig/sum(pco$CA$eig)*100

  df<-data.frame(x=summary(pco)$sites[,axis1],y=summary(pco)$sites[,axis2],Group=map[,color.group],ID=map$ID)
  col1<-figure.colors[factor(map[,color.group])]

  #t.test across MDS1 and MDS2
  df.paired<- df %>% group_by(ID) %>% filter (n()==2) %>% as.data.frame()
  results <- data.frame(Group = c('HC-T1','HC-T2','T1-T2'))

  for (variable in c('x','y')){
    compare.groups1<-t.test(df[,variable][df$Group=="HC"],
                            df[,variable][df$Group=="T1"])
    compare.groups2<-t.test(df[,variable][df$Group=="HC"],
                            df[,variable][df$Group=="T2"])
    compare.groups3<-t.test(df.paired[,variable][df.paired$Group=="T1"],
                            df.paired[,variable][df.paired$Group=="T2"],paired = TRUE)
    results[,variable]<-p.adjust(c(compare.groups1$p.value,compare.groups2$p.value,compare.groups3$p.value),method = 'BH')
  }

  pdf(pdf.name,width = 6,height=6)

  #boxplot for MDS2 (y axis)
  par(fig=c(0.04,0.33,0.2,1),new=FALSE)
  boxplot(df$y~df$Group,col="white",ann = FALSE,boxwex=0.8,outline=FALSE,yaxt="n",xaxt="n",frame.plot=FALSE)
  points(df$y~jitter(as.numeric(as.factor(df$Group)),0.4),col=col1,cex=0.3,pch=16)

  #PCOA
  par(fig=c(0.27,1,0.23,1),new=TRUE)
  pcoa=ordiplot(pco,choices=c(axis1,axis2),display="sites",type="none",cex.lab=0.8,
                xlab=paste0("MDS",axis1, " (",format(percent[axis1],digits = 3),"%)"),
                ylab=paste0("MDS",axis2," (",format(percent[axis2],digits = 3),"%)"),
                cex.main=1,main=plot.title)
  points(pcoa,"sites",col=adjustcolor(col1,alpha.f = 0.5),pch=16,cex=0.7)

  if(show.ellipses)
    ordiellipse(pcoa, factor(map[,color.group]), kind="se", conf=0.95, lwd=1.5,
                draw = "lines", col=cols,
                show.groups=levels(factor(map[,color.group])),
                label=T,font=2,cex=0.6)

  if(show.legend)
    legend("bottomleft",legend = legend.names,
           col=figure.colors,pch=16,cex=0.7)

  #Boxplot for MDS1, x axis
  par(fig=c(0.25,1,0,0.39),new=TRUE)
  boxplot(df$x~df$Group,col="white",ann = FALSE,boxwex=0.7,outline=FALSE,yaxt="n",,horizontal =TRUE,xaxt="n",frame.plot=FALSE)
  points(jitter(as.numeric(as.factor(df$Group)),0.4) ~ df$x,col=col1,cex=0.3,pch=16)

  dev.off()
  return(results)
}

#Boxplots and t.test for phenotypes

phenotype.boxplot<-function(map,variable,y.title,plot.title,show.legend=FALSE){


  map.paired<-map %>% group_by(ID,) %>% filter(n()==2) %>% as.data.frame()

  if(length(unique(map$Location))==1 && unique(map$Location)=="Denver" && !variable %in% c('BMI','shannon_taxonomy','shannon_pathways')){

    compare.groups<-t.test(map.paired[,variable][map.paired$Type=="T1"],
                           map.paired[,variable][map.paired$Type=="T2"],
                           paired=TRUE)
    pvals = paste0('T1 vs.T2 p = ',format(compare.groups$p.value,digits = 3))
    data=map.paired

  } else {

    compare.groups1<-t.test(map[,variable][map$Type=="HC"],
                            map[,variable][map$Type=="T1"])
    compare.groups2<-t.test(map[,variable][map$Type=="HC"],
                            map[,variable][map$Type=="T2"])
    compare.groups3<-t.test(map.paired[,variable][map.paired$Type=="T1"],
                            map.paired[,variable][map.paired$Type=="T2"],paired = TRUE)

    pvals_adjusted <-p.adjust(c(compare.groups1$p.value,compare.groups2$p.value,compare.groups3$p.value),method = 'BH')

    pvals<-paste0('non-ED vs. T1 adjusted p = ',format(pvals_adjusted[1],digits = 3),
                  '\nnon-ED vs. T2 adjusted p = ',format(pvals_adjusted[2],digits = 3),
                  '\nT1 vs. T2 adjusted p = ',format(pvals_adjusted[3],digits = 3))
    data=map
  }
  cols<-brewer.pal(3, "Dark2")
  plot<-ggplot(data,aes(x=Type,y=data[,variable],color=Type))+
    geom_boxplot(outlier.shape=NA)+
    geom_jitter(position = position_jitter(width = 0.1),size=0.8)+
    scale_color_manual(values=cols)+scale_x_discrete(labels = names.type)+
    labs(x="",y=y.title,title = plot.title,
         subtitle =pvals )+
    theme_classic(10)
  if(!show.legend)
    plot<-plot+theme(legend.position = "none")
  return(plot)
}

#Boxplot for all taxa comparing HC, T1, T2
get.box.plots<-function(map,otu,result.test,FDR=0.1,order.by.column,is.taxonomy){

  data.all<-cbind(map,otu)
  sig.result<-result.test %>% filter (.data[[order.by.column]] < FDR )

  if(length(sig.result[,order.by.column])==0){
    cat("Boxplots comparing groups: No Significant hit for",order.by.column,"at FDR",FDR,"\n")
  }else{

    sig.result<-sig.result %>% arrange(.data[[order.by.column]])
    plots<-apply(sig.result,1,function(x) box.plot(data.all,x,show.legend=FALSE,is.taxonomy))
    return(plots)
  }
}

#Boxplot for each taxon comparing HC, T1, T2
box.plot<-function(data.all,pval.taxa,show.legend=TRUE,is.taxonomy){

  if(is.taxonomy){
    y.title = expression('log'[10]~'Normalized Count')
  } else{
    y.title ='Normalized Count'
  }

  taxaName<-as.character(pval.taxa[1])
  pval.HC.T1<-format(as.numeric(pval.taxa[4]),digit=2)
  pval.HC.T2<-format(as.numeric(pval.taxa[7]),digit=2)
  pval.T1.T2<-format(as.numeric(pval.taxa[10]),digit=2)

  p<-ggplot(data=data.all,aes(y=data.all[,taxaName],x=Type))+
    geom_boxplot(outlier.shape=NA)+theme_classic()+
    labs(title=taxaName,y=y.title,x="",
         subtitle = paste("non-ED vs. T1 adjusted p =",pval.HC.T1,"\nnon-ED vs. T2 adjusted p =",pval.HC.T2,"\nT1 vs. T2 adjusted p =",pval.T1.T2))+
    geom_jitter(position = position_jitter(width = 0.1),aes(col=Type),size=0.8,alpha=0.8 )+
    scale_color_brewer(palette = "Dark2")+
    scale_x_discrete(labels = names.type)

  if(!show.legend){
    p<-p+theme(legend.position = "none")
  }
  return(p)
}


scatter.plot<-function(map,otu,taxa,variable,xlab,pvalue,legend.show=FALSE){

  data.all<-cbind(map,otu)

  if(length(unique(data.all$Location))>1){
    colors<-c("#1F78B4","#F0027F")
  } else if (length(unique(data.all$Location))==1 && unique(data.all$Location)=="Denver"){
    colors<-"#1F78B4"
  }else {
    colors<-"#F0027F"
  }

  p<-ggplot(data=data.all,aes(x=data.all[,variable],y=data.all[,taxa]))+
    geom_point(aes(color=Location),alpha=0.8)+
    labs(x=xlab,y=expression('log'[10]~'normalized count'),
         title=taxa,subtitle = paste('adjusted p =',format(pvalue,digits = 3)))+
    scale_color_manual(values =colors)+
    theme_classic()

  if(!legend.show){
    p<-p+theme(legend.position = "none")
  }
  p
}

#Scatter plots for all variables
get.scatter.plots<-function(map,otu,variable,variable.name,result.test,FDR=0.1,legend.show){

  if(!is.data.frame(result.test)){
    result.test<-as.data.frame(result.test)
  }
  sig.result<-result.test %>% filter (result.test$`adjusted.p-variable` < FDR)

  if(length(sig.result$`adjusted.p-variable`)==0){
    cat("Scatter plot for variable.name: No significant hit for",variable.name," at FDR",FDR,"\n")
  }else{
    sig.result<-sig.result %>% arrange(sig.result$`adjusted.p-variable`)
    p.values<-sig.result$`adjusted.p-variable`
    taxa=rownames(sig.result)

    plots<-lapply(taxa, function(x) scatter.plot(map=map,otu=otu,
                                                 taxa=x,variable=variable,xlab=variable.name,
                                                 pvalue=p.values[which(taxa==x)],legend.show = legend.show))

    return(plots)
  }
}

#Create log10 pvalue versus log10 p-value
plot.pvals<-function(pval1,pval2,name1,name2,plot.title){

  df<-data.frame(pval1,pval2)
  df_sig <- df[(df$pval1 < log10(0.05) | pval1 > -log10(0.05)) &
                 (df$pval2 < log10(0.05) | pval2 > -log10(0.05)),]

  spearman.rho<-correlation.pvals(df)[[1]]
  spearman.p<-correlation.pvals(df)[[2]]

  plot<-ggplot(data=df,aes(x=pval1,y=pval2))+
    geom_point(size=0.5,color=alpha("black",0.5))+
    geom_hline(yintercept = 0,linetype="dashed", color = "red")+
    geom_vline(xintercept = 0,linetype="dashed", color = "red")+
    labs(x=name1,y=name2,title=plot.title,
         subtitle = paste0("Spearman coefficient = ",format(spearman.rho,digits =3),
                           "\nadj p = ",format(spearman.p,digits = 3)))+
    theme_classic()+
    geom_point(data=df_sig, aes(color = "red"),size=0.5)+
    theme(legend.position = 'none')

  plot
}

#Create a bar plot showing the effect size from adonis test
plot.adonis<-function(adonis.result,variable.names,file.name,FDR=0.1,show.legend=TRUE){

  adonis.result<-as.data.frame(adonis.result)
  adonis.result$variables<-rownames(adonis.result)
  adonis.result$sig<-ifelse(adonis.result$adjusted.pvalue<FDR,"significant","non-significant")
  adonis.result$names<-variable.names
  plot<-ggplot(data=adonis.result,aes(x=R2,y=variable.names))+geom_bar(stat="identity",aes(fill=sig),width = 0.5)+
    scale_fill_manual(values = c("gray","red"))+labs(y="",x=expression('R'^2))+theme_classic()
  if(!show.legend)
    plot<-plot+theme(legend.position = "none")
  pdf(file.name,width = 6,height = 5)
  print(plot)
  dev.off()
  plot
}

