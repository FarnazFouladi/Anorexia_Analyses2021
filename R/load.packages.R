packages<-c('circlize','ComplexHeatmap','vegan','ggplot2',
            'RColorBrewer','cowplot','nlme','dplyr','Hmisc','gridExtra',
            'ggpubr','rstatix')

installed_packages <- packages %in% rownames(installed.packages())

if(any(installed_packages) == FALSE){
  install.packages(packages[!installed_packages])
}

#Load libraries
invisible(lapply(packages, library, character.only = TRUE))
