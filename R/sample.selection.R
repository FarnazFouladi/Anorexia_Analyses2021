select.samples<-function(map,otu,site=NULL,disease=NULL){
  
  if(!is.null(site) & !is.null(disease)){
    map_sub<-map[map$Location==site & map$Type %in% disease,]
    otu_sub<-otu[map$Location==site & map$Type %in% disease,]
  }else if (!is.null(site) & is.null(disease)){
    map_sub<-map[map$Location==site,]
    otu_sub<-otu[map$Location==site,]
  }else if(is.null(site) & !is.null(disease)){
    map_sub<-map[map$Type %in% disease,]
    otu_sub<-otu[map$Type %in% disease,]
  }else{
    map_sub<-map
    otu_sub<-otu
  }
  list(map_sub,otu_sub)
}


get.paired.samples<-function(map,otu){
  pairedID<-intersect(map$ID[map$Type=="T1"],map$ID[map$Type=="T2"])
  map_paired<-map %>% filter(ID %in% pairedID) %>% group_by(Type) %>% arrange(ID,.by_group = TRUE)
  otu_paired<-otu[map_paired$SampleID,]
  list(map_paired,otu_paired)
}

