t = "Pathway"
t.test.result<-read.table(paste0("output/",t,"/",t,"_t-test_All.txt"),sep = "\t",header = TRUE,row.names = 1,check.names = FALSE,quote = "")
result <- t.test.result
result <- result %>% filter(adj.p.vals.HC.T1 < 0.05 | adj.p.vals.HC.T2 < 0.05 | adj.p.vals.T1.T2 < 0.05)
result$adj.p.vals.HC.T1[result$adj.p.vals.HC.T1 > 0.05]  <- 1
result$adj.p.vals.HC.T2[result$adj.p.vals.HC.T2 > 0.05]  <- 1
result$adj.p.vals.T1.T2[result$adj.p.vals.T1.T2 > 0.05]  <- 1


result$logHCT1 <- getlog10p(result$adj.p.vals.HC.T1,result$t.HC.T1)
result$logHCT2 <- getlog10p(result$adj.p.vals.HC.T2,result$t.HC.T2)
result$logT1T2 <- getlog10p(result$adj.p.vals.T1.T2,result$t.T1.T2)


annotated_path<- readxl::read_xlsx('output/Pathway/significant_pathways.xlsx',na = 'NA')
annotated_path$Pathways_new <- as.character(sapply(annotated_path$Pathways, function(x) strsplit(x,':')[[1]][1]))

annotated_path$level1 <- sapply(as.character(annotated_path$Class),function(x) {if (x!='Superpathways') strsplit(x,';')[[1]][1] else return('Superpathways')})
annotated_path$level2 <- sapply(as.character(annotated_path$Class),function(x) {if (x!='Superpathways') strsplit(x,';')[[1]][2] else return('Superpathways')})
annotated_path$Pathways_new <- paste0('L1__Pathways;','L2__',annotated_path$level1,';', 'L3__',annotated_path$level2,';','L4__',annotated_path$Pathways_new)

merged_df <- merge(result, annotated_path, by.x = 'taxa', by.y = 'Pathways', all.x = TRUE)

merged_df <- merged_df %>% mutate(pathways = Pathways_new) %>%
  filter(!is.na(pathways))

library(metacoder)
obj <- parse_tax_data(merged_df,
                      class_cols = "pathways", # the column that contains taxonomic information
                      class_sep = ";", # The character used to separate taxa in the classification
                      class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
                      class_key = c(tax_rank = "info", # A key describing each regex capture group
                                    tax_name = "taxon_name"))

taxon_id <- names(obj$taxon_names())
taxon_name <- as.character(obj$taxon_names())

df <- data.frame(taxon_id, taxon_name)
merged_df$pathways <- sapply(merged_df$pathways, function(x) strsplit(x,'L4__')[[1]][2])
df1 <- merge(df, merged_df, by.x = 'taxon_name',by.y ='pathways',all.x = TRUE )
df1 <- df1[match(taxon_id,df1$taxon_id),]
df1$logHCT1[is.na(df1$logHCT1)] <- 0
df1$logHCT2[is.na(df1$logHCT2)] <- 0
df1$logT1T2[is.na(df1$logT1T2)] <- 0
df1 <- as_tibble(df1)
obj$data$tax_pvals <- df1


pdf('Figures/Figure2/tree_all.pdf',20,20)

set.seed(2)
plot1 <- heat_tree(obj,
                  node_label = obj$taxon_names(),
                  node_color = obj$data$tax_pvals$logHCT1,
                  node_color_interval = c(-15, 15),
                  node_color_range = c("blue", "gray", "red"),
                  node_color_axis_label = 'log10 p-value',
                  node_size_axis_label = NULL,
                  node_label_size = 0.005,
                  title  = 'HC-T1')
print(plot1)


plot2 <- heat_tree(obj,
                   node_label = obj$taxon_names(),
                   node_color = obj$data$tax_pvals$logHCT2,
                   node_color_interval = c(-15, 15),
                   node_color_range = c("blue", "gray", "red"),
                   node_color_axis_label = 'log10 p-value',
                   node_size_axis_label = NULL,
                   node_label_size =0.005,
                   title  = 'HC-T2')

print(plot2)
plot3 <- heat_tree(obj,
                   node_label = obj$taxon_names(),
                   node_color = obj$data$tax_pvals$logT1T2,
                   node_color_interval = c(-15, 15),
                   node_color_range = c("blue", "gray", "red"),
                   node_color_axis_label = 'log10 p-value',
                   node_size_axis_label = NULL,
                   node_label_size =0.005,
                   title  = 'T1-T2')
print(plot3)
dev.off()



