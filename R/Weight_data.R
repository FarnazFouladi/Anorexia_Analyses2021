meta <- read.table('data/meta_master_noweight.txt',sep='\t',header=TRUE)
weight <- read.table('data/weight_data.txt',sep='\t',header=TRUE)
weight$SampleID <- sapply(weight$SampleID,function(x) strsplit(x,' ')[[1]][1])
weight <- weight[,c(1,2)]
weight$ID <- sapply(weight$SampleID,function(x) substr(x,1,6))
weight$Type <-  sapply(weight$SampleID,function(x) substr(x,7,8))

df <- weight %>% filter(Type!='HC') %>% group_by(ID) %>% filter (n()==2) %>% as.data.frame()
changeInWeight <- sapply(seq(1,nrow(df),2),function(x) df$Weight.kg.[x+1] - df$Weight.kg.[x])
df$changeInWeight  <- rep(changeInWeight ,each = 2)
df <- df %>% select(c("SampleID","changeInWeight"))

weight_withChange <- merge(weight,df,by='SampleID',all=TRUE)
weight_withChange <- weight_withChange %>% select(-c('Type',ID))

meta_merged <- merge(meta,weight_withChange,by='SampleID',all = TRUE)
meta_merged$changeInWeightPerDay <- meta_merged$changeInWeight/meta_merged$TimeBetweenRecovery

write.table(meta_merged,'data/meta_master.txt',sep='\t',quote = FALSE)



