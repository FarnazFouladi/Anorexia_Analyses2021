# If paired test is needed pass the results of getOnlyPairedSamples
# t.test (paired between T1 and T2, unpaired between HC and T1,unpaired between HC and T2)
compare.Groups<-function(map,otu,group1,group2,paired=FALSE){

  if(paired){
    map<-get.paired.samples(map,otu)[[1]]
    otu<-get.paired.samples(map,otu)[[2]]
  }
  p.vals<-apply(otu,2,function(x) t.test(x[map$Type==group1],x[map$Type==group2],paired = paired)$p.value)
  t.stat<-apply(otu,2,function(x) t.test(x[map$Type==group1],x[map$Type==group2],paired = paired)$statistic)

  adjusted.p<-p.adjust(p.vals,method = "BH")

  df<-data.frame(taxa=names(p.vals),t.stat,p.vals,adjusted.p)
  colnames(df)=c("taxa",paste0("t.",group1,".",group2)
                 ,paste0("p.vals.",group1,".",group2),
                 paste0("adj.p.vals.",group1,".",group2))
  return(df)
}

#Adonis test
perform.adnois<-function(map,otu,variable){
  df<-cbind(otu,variable=map[,variable])
  df1<-na.omit(df)
  fit<-adonis(df1[,1:ncol(otu)] ~ df1[,"variable"])
  stats<-c(fit$aov.tab[1,5],fit$aov.tab[1,6])
  names(stats)<-c("R2","p-value")
  stats
}

#Adonis test on all variables
perform.adonis.all.vars<-function(otu,map,variables,file.Name){
  p.vals<-lapply(variables,perform.adnois,otu=otu,map=map)
  p.vals.matrix<-matrix(unlist(p.vals),nrow=length(variables),ncol=2,byrow = TRUE)
  p.vals.matrix<-cbind( p.vals.matrix,p.adjust(p.vals.matrix[,2],method = "BH"))
  rownames(p.vals.matrix)<-variables
  colnames(p.vals.matrix)<-c("R2","pvalue","adjusted.pvalue")
  write.table(p.vals.matrix,file.Name,quote = FALSE,sep="\t")
  return(p.vals.matrix)
}


#Transform p-values -log10(pvalue) * sign of coefficient
#If a coefficient is positive, the log10 p-value will be positive.
getlog10p<-function(pval,coefficient){

  log10p<-sapply(1:length(pval), function(x){

    if(!is.na(pval[x])){
      if (coefficient[x]>0) return(-log10(pval[x]))
      else return(log10(pval[x]))
    }else{
      return(NA)
    }
  })
  log10p
}

#Correlation between p-values
correlation.pvals<-function(df){

  myList<-list()
  test<-cor.test(df[,"pval1"],df[,"pval2"],method = "spearman")
  myList[[1]]<-test$estimate
  myList[[2]]<-test$p.value
  return(myList)
}



