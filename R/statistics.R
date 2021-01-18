require(nlme)
library(dplyr)

# If paired test is needed pass the results of getOnlyPairedSamples
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

#Mixed linear model (Add t statitics and figure out a way for within cohort)

MLM<-function(taxa,variable,treatment,cohort,ID){

  myData<-data.frame(bug=taxa,variable =variable,
                     treatment=factor(treatment),cohort=factor(cohort),ID)

  fit<-tryCatch({summary(lme(bug ~ variable * treatment , random = ~ 1 | cohort/ID, data=myData,na.action = na.omit))},
                error=function(e){cat("ERROR :",conditionMessage(e), "\n")
                  return(NA)})

  if (!is.na(fit)){
    pvals<-c(fit$tTable[2,5],fit$tTable[3,5],fit$tTable[4,5],fit$tTable[2,1])
  }else{
    pvals<-c(NA,NA,NA,NA)
  }

  names(pvals)<-c("p-variable","p-treatment","p-interaction","slope")
  pvals
}

perform.MLM<-function(otu,map,variable){

  p.vals<-apply(otu,2,function(i) MLM(taxa=i,variable =map[,variable],treatment = map[,"Type"],ID=map[,"ID"],cohort = map[,"Location"] ))
  p.vals<-t(p.vals)
  adjusted.p<-apply(p.vals[,1:3],2,p.adjust,method="BH")
  colnames(adjusted.p)<-paste0("adjusted.",colnames(p.vals[,1:3]))
  p.vals<-cbind(p.vals,adjusted.p)
  return(p.vals)
}

perform.MLM.all.vars<-function(map,otu,variables,file.Name){
  p.vals<-lapply(variables,perform.MLM,otu=otu,map=map)
  mapply(write.table,x=p.vals,file=file.Name,quote = FALSE,sep="\t")
  names(p.vals)<-variables
  p.vals
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

perform.adonis.all.vars<-function(otu,map,variables,file.Name){
  p.vals<-lapply(variables,perform.adnois,otu=otu,map=map)
  p.vals.matrix<-matrix(unlist(p.vals),nrow=length(variables),ncol=2,byrow = TRUE)
  p.vals.matrix<-cbind( p.vals.matrix,p.adjust(p.vals.matrix[,2],method = "BH"))
  rownames(p.vals.matrix)<-variables
  colnames(p.vals.matrix)<-c("R2","pvalue","adjusted.pvalue")
  write.table(p.vals.matrix,file.Name,quote = FALSE,sep="\t")
  return(p.vals.matrix)
}


#Transform p-values
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



