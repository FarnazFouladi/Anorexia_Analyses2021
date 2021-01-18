library(vegan)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

pco.plot<-function(map,otu,color,shape,plot.title,axis1=1,axis2=2,
                   show.legend=TRUE,show.stat=TRUE){

  fit<-adonis(otu ~ factor(map[,color]),method = "bray")
  pval=fit$aov.tab$`Pr(>F)`[1]
  R2=fit$aov.tab$R2[1]
  pval_round=format(pval,digits = 3)
  R2_round=format(R2,digits = 3)

  pco<-capscale(otu ~ 1, distance = "bray")
  percent<-pco$CA$eig/sum(pco$CA$eig)*100

  if(is.null(shape))
    df<-data.frame(x=pco$CA$u[,axis1],y=pco$CA$u[,axis2],Group=map[,color])
  else
    df<-data.frame(x=pco$CA$u[,axis1],y=pco$CA$u[,axis2],Group=map[,color],Site=map[,shape])

  if(length(unique(df$Site))==1 && unique(df$Site)=="UNC"){
    point.shape=1
  } else if (length(unique(df$Site))==1 && unique(df$Site)=="Denver"){
    point.shape=16
  } else{
    point.shape=c(16,2,1)
  }

  p<-ggplot(data=df, aes(x, y))

  if(is.null(shape))
    p<-p+geom_point(aes(colour=Group),alpha=0.8)+scale_color_brewer(palette = "Dark2")
  else
    p<-p+geom_point(aes(colour=Group,shape=Site),alpha=0.8)+
    scale_color_brewer(palette = "Dark2")+
    scale_shape_manual(values = point.shape)

  p<-p+xlab(paste0("PC",axis1," (",format(percent[axis1],digit=3),"%)")) +
    ylab(paste0("PC",axis2," (",format(percent[axis2],digit=3),"%)")) +
    ggtitle(plot.title) +
    theme_minimal()+
    stat_ellipse(aes(colour=Group), show.legend=F, type="t", level=.6)

  if (!show.legend){
    p<-p+theme(legend.position = "none")
  }
  if(show.stat){
    p<-p+labs(subtitle =bquote('p =' ~ .(pval_round)~'R'^2~'='~.(R2_round)))
  }
  return(p)
}



get.box.plots<-function(map,otu,result.test,FDR=0.1,order.by.column){

  data.all<-cbind(map,otu)
  sig.result<-result.test %>% filter (.data[[order.by.column]] < FDR )

  if(length(sig.result[,order.by.column])==0){
    cat("No Significant hit for",order.by.column,"at FDR",FDR)
  }else{

    sig.result<-sig.result %>% arrange(.data[[order.by.column]])
    plots<-apply(sig.result,1,function(x) box.plot(data.all,x,show.legend=FALSE))
    return(plots)
  }
}


box.plot<-function(data.all,pval.taxa,show.legend=TRUE){

  taxaName<-as.character(pval.taxa[1])
  pval.HC.T1<-format(as.numeric(pval.taxa[4]),digit=2)
  pval.HC.T2<-format(as.numeric(pval.taxa[7]),digit=2)
  pval.T1.T2<-format(as.numeric(pval.taxa[10]),digit=2)

  p<-ggplot(data=data.all,aes(y=data.all[,taxaName],x=Type))+
    geom_boxplot()+theme_bw()+
    labs(title=taxaName,y=expression('log'[10]~'Normalized Count',),x="",
         subtitle = paste("HC vs. T1 p =",pval.HC.T1,"\nHC vs. T2 p =",pval.HC.T2,"\nT1 vs. T2 p =",pval.T1.T2))+
    geom_jitter(position = position_jitter(width = 0.1),aes(col=Type),size=0.8,alpha=0.8 )+
    scale_color_brewer(palette = "Dark2")

  if(!show.legend){
    p<-p+theme(legend.position = "none")
  }
  return(p)
}


scatter.plot<-function(map,otu,taxa,variable,xlab,pvalue,legend.show=FALSE){

  data.all<-cbind(map,otu)

  p<-ggplot(data=data.all,aes(x=data.all[,variable],y=data.all[,taxa]))+
    geom_point(aes(color=Location),alpha=0.8)+
    labs(x=xlab,y=expression('log'[10]~'normalized count'),
         title=taxa,subtitle = paste('p =',format(pvalue,digits = 3)))+
    scale_color_manual(values = c("green","hotpink","purple"))+
    theme_bw()

  if(!legend.show){
    p<-p+theme(legend.position = "none")
  }
  p
}


get.scatter.plots<-function(map,otu,variable,variable.name,result.test,FDR=0.1,legend.show){

  if(!is.data.frame(result.test)){
    result.test<-as.data.frame(result.test)
  }
  sig.result<-result.test %>% filter (result.test$`adjusted.p-variable` < FDR)

  if(length(sig.result$`adjusted.p-variable`)==0){
    cat("No significant hit at FDR",FDR)
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


plot.pvals<-function(pval1,pval2,name1,name2,plot.title){

  df<-data.frame(pval1,pval2)
  spearman.rho<-correlation.pvals(df)[[1]]
  spearman.p<-correlation.pvals(df)[[2]]

  plot<-ggplot(data=df,aes(x=pval1,y=pval2))+geom_point(size=0.5,color=alpha("black",0.5))+
    geom_hline(yintercept = 0,linetype="dashed", color = "red")+
    geom_vline(xintercept = 0,linetype="dashed", color = "red")+
    labs(x=name1,y=name2,title=plot.title,
         subtitle = paste0("Spearman coefficient = ",format(spearman.rho,digits =3),
                           "\nadj p = ",format(spearman.p,digits = 3)))+
    theme_classic()
  plot
}

plot.adonis<-function(adonis.result,variable.names,file.name,FDR=0.1,show.legend=TRUE){

  adonis.result<-as.data.frame(adonis.result)
  adonis.result$variables<-rownames(adonis.result)
  adonis.result$sig<-ifelse(adonis.result$adjusted.pvalue<FDR,"significant","non-significant")
  adonis.result$names<-variable.names
  plot<-ggplot(data=adonis.result,aes(x=R2,y=variable.names))+geom_bar(stat="identity",aes(fill=sig),width = 0.5)+
    scale_fill_manual(values = c("gray","red"))+labs(y="",x=expression('R'^2))+theme_bw()
  if(!show.legend)
  plot<-plot+theme(legend.position = "none")
  pdf(file.name,width = 6,height = 5)
  print(plot)
  dev.off()
  plot
}

