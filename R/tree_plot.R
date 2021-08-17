library(ggplot2)
library(ggpubr)
library(ggtree)
library(RColorBrewer)

path1_1=read.csv(file="significant_pathways.csv",header=T)
path1=path1_1
path1[,1]=sapply(strsplit(path1_1[,1],": "),"[[",1)
path2=paste(path1[,2],path1_1,sep=";")

a=matrix(nrow=nrow(path1),ncol=7)
for (n in 1:nrow(path1)){
  a1=strsplit(path1[n,2],";")[[1]]
  a[n,]=rep("none",7)
  a[n,1:length(a1)]=a1
  a[n,7]=path1[n,1]
}
a=gsub(";","--",a)
a=gsub(":","-",a)
a=gsub(" ","",a)
a=gsub("/","",a)
a=gsub("\\'","",a)
a=gsub(",","",a)
a=gsub("\\s","",a)


b1=apply(a,1,function(i){paste(i,collapse="--")})
b2=unique(apply(a,1,function(i){paste(i[1:6],collapse="--")}))
b3=unique(apply(a,1,function(i){paste(i[1:5],collapse="--")}))
b4=unique(apply(a,1,function(i){paste(i[1:4],collapse="--")}))
b5=unique(apply(a,1,function(i){paste(i[1:3],collapse="--")}))
b6=unique(apply(a,1,function(i){paste(i[1:2],collapse="--")}))
b7=unique(apply(a,1,function(i){paste(i[1],collapse="--")}))
bn1=c(b1,b2,b3,b4,b5,b6,b7)


ln=sapply(strsplit(bn1,"--"),length)
names1=bn1
names7=bn1[ln==7]
l1=list()
for (n in 7:2){
  l2=list()
  if(n==7){
    for (m in names1[which(ln==n-1)]){
      names2=names1[ln==n]
      l1[m]=paste0("(",paste(names2[grep(m,names2)],collapse=","),")",m)
    }
  }else{
    for (m in names1[which(ln==n-1)]){
      names2=names1[ln==n]
      l2[m]=paste0("(",paste(l1[names2[grep(m,names2)]],collapse=","),")",m)
    }
    l1=l2
  }
}
l2=paste0("(",paste(l1,collapse=","),")Pathway;")
tree <- ape::read.tree(text = l2)
p=ggtree::ggtree(tree, layout='circular',edge.length=NULL)

node_num=vector()
tip_num=vector()
for (i in 1:length(b6)){
  print(i)
  nodem=try(ggtree::MRCA(tree,tree$tip.label[grep(as.character(b6[i]),tree$tip.label)]))
  if(class(nodem)=="try-error"){
    node_num[i]=NA
    tip_num[i]=NA
  }else if (length(grep(as.character(b6[i]),tree$tip.label))<4){
    node_num[i]=NA
    tip_num[i]=length(grep(as.character(b6[i]),tree$tip.label))
  }else{
    node_num[i]=nodem
    tip_num[i]=length(grep(as.character(b6[i]),tree$tip.label))
  }
}

dd3=data.frame(cbind(b6,node_num,tip_num))
dd3$node_num=as.numeric(dd3$node_num)
dd3$name1=sapply(strsplit(dd3$b6,"--"),"[[",2)
dd3=data.frame(dd3[dd3$name1!="none" & !is.na(dd3$node_num),])
dd3=dd3[which(!is.na(dd3$node_num)),]
dd3=dd3[order(dd3$node_num),]

dd3$name1
dd3$Category=c("Amine and Polyamine\nBiosynthesis", "Amino Acid\nBiosynthesis" ,"Carbohydrate\nBiosynthesis",
               "Cell Structure\nBiosynthesis" ,"Cofactor Carrier and\nVitamin Biosynthesis", "Fatty Acid and\nLipid Biosynthesis",
               "Nucleoside and Nucleotide\nBiosynthesis","Secondary Metabolite\nBiosynthesis" , "Amine and Polyamine\nDegradation",
               "Carbohydrate\nDegradation" , "Carboxylate\nDegradation" ,"Nucleoside and Nucleotide\nDegradation"  ,
               "Fermentation" ,"Glycolysis","TCAcycle" )

p=p + geom_hilight(data=dd3,mapping=aes(fill=Category, node=node_num),alpha=0.2,extend=0.5) +
  scale_fill_manual(values=rainbow(nrow(dd3)))

txt1=c(1.5,1,1.2,1,1,2.5,1,2.8,2.5,1,2.5,1,1,1,1.5)
p1=p
for(n in 1:nrow(dd3)){
  p1=p1 + ggtree::geom_cladelabel(node = dd3$node_num[n], label = dd3$Category[n],barsize =0,angle = 'auto',horizontal = F,offset.text = txt1[n],hjust = "center", fontsize = 3)
}


ggt=list()
compares=list()
compares[[1]]=c("nonED","T1")
compares[[2]]=c("nonED","T2")
compares[[3]]=c("T1","T2")

title1=list()
title1[[1]]="All"
title1[[2]]="CEED"
title1[[3]]="ACUTE"

re=read.csv(file="results_pathway.csv",header=T,row.names=1)
for (i in 1:9){
  re1=re[which(re[,i*3]<0.1),(i*3-2):(i*3)]
  re1$name=b1[match(rownames(re1),path1_1[,1])]
  re1$group1=-sign(re1[,1])
  m=i%%3
  if(m==0){
    m=3
  }
  re1$group=compares[[m]][factor(re1$group1)]
  col1n=c("red","blue","goldenrod1")[m-4]

  re1=re1[,c(4,1:3,5:6)]
  main1=paste(title1[[ceiling(i/3)]],paste(compares[[m]],collapse = " vs "))
  ggt[[i]]=p1 %<+% re1 + geom_point(aes(color=group), alpha=1)+ scale_colour_manual(values=col1n,na.translate = F)+ggtitle(main1)
}
ggexport(plotlist = ggt, filename = "tree1.pdf",width=10,height=10)
