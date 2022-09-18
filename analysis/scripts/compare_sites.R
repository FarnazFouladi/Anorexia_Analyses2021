# Differences between in microbiome between the two sites

# Perform t-test

source("../resources/PrepareDataForFigures.R")

for ( t in c("genus","species","pathway")){

  # Create output directory
  output_dir <- paste0('../output/compare_sites/',capitalize(t))
  dir.create(output_dir,recursive = T,showWarnings = F)

  # get the normalized table
  df_temp <- get(t)

  # Check the order of count table and map
  stopifnot(all.equal(rownames(df_temp),map$SampleID))

  for (type in unique(map$Type)){

    # Subset data to the disease type
    map_sub <- map %>% filter(Type == type)
    df_sub <- df_temp[map_sub$SampleID,]
    df_sub <- df_sub[,colSums(df_sub)>0]
    stopifnot(all.equal(map_sub$SampleID,rownames(df_sub)))

    results <- matrix(NA, ncol = 3, nrow = ncol(df_sub), dimnames = list(colnames(df_sub),c('t_statistics','p_value','adjusted_p_value')))

    for (taxa in colnames(df_sub)){

      # t-test
      fit <- t.test(df_sub[,taxa] ~ map_sub[,'Location'])
      results[taxa,1]  <- fit$statistic
      results[taxa,2] <- fit$p.value
    }

    # adjusted p-value
    results[,3] <- p.adjust(results[,2],method = "BH")

    # Create feature column, order by adjusted p-value
    results <- results %>% as.data.frame() %>%
      tibble::rownames_to_column('Feature') %>%
      arrange(adjusted_p_value)

    # Adnonis test
    adonis.result <- perform.adonis.all.vars(otu = df_sub,map = map_sub,variables = 'Location',file.path(output_dir,paste0(type,'_',t,'_adonis.txt')))
    write.table(results, file.path(output_dir,paste0(type,'_',t,'_compare_sites.txt')),sep = '\t',row.names = F)

  }
}

# The only significant difference is Bifidobacterium animalis between Denver T2 and UNC
species_name = 'Bifidobacterium animalis'

df_merged <- cbind(species = species[,'Bifidobacterium animalis'],map)
plot_bifidobacterium <- ggplot(df_merged, aes(y = species, x = Location))+
  geom_boxplot()





