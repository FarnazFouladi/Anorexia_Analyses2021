# Figure1 -----------------------------------------------------------------

dir.name <- 'Figures/Figure1'
if(!dir.exists(dir.name)){
  dir.create(dir.name,recursive = TRUE,showWarnings = FALSE)
}

source('R/PrepareDataForFigures.R')

##########Taxonomy Ordination Plot##########

#CEED and ACUTE

p1<-pco.plot(map,species,color.group = "Type",
             plot.title = "Taxonomy-CEED and ACUTE",figure.colors = cols.type,
             show.legend = TRUE,legend.names = names.type,
             pdf.name = "Figures/Figure1/Species_PCO1&2_all.pdf")

p2<-pco.plot(map,species,color.group = "Type",axis1=3,axis2 = 4,
             plot.title = "Taxonomy-CEED and ACUTE",figure.colors = cols.type,
             legend.names = names.type,
             pdf.name = "Figures/Figure1/Species_PCO3&4_all.pdf")

#CEED
p3<-pco.plot(map.unc,species.unc,color.group = "Type",
             plot.title = "Taxonomy-CEED",figure.colors = cols.type,
             show.legend = TRUE,legend.names = names.type,
             pdf.name = "Figures/Figure1/Species_PCO1&2_CEED.pdf")

p4<-pco.plot(map.unc,species.unc,color.group = "Type", axis1=3,axis2 = 4,
             plot.title = "Taxonomy-CEED",figure.colors = cols.type,
             show.legend = TRUE,legend.names = names.type,
             pdf.name = "Figures/Figure1/Species_PCO3&4_CEED.pdf")

#ACUTE
p5<-pco.plot(map.denver,species.denver,color.group = "Type",
             plot.title = "Taxonomy-ACUTE",figure.colors = cols.type,
             show.legend = TRUE,legend.names = names.type,
             pdf.name = "Figures/Figure1/Species_PCO1&2_ACUTE.pdf")

p6<-pco.plot(map.denver,species.denver,axis1=3,axis2 = 4,color.group = "Type",
             plot.title = "Taxonomy-ACUTE",figure.colors = cols.type,
             show.legend = TRUE,legend.names = names.type,
             pdf.name = "Figures/Figure1/Species_PCO3&4_ACUTE.pdf")


##########Taxonomy Ordination Plot##########
#CEED and ACUTE
p7<-pco.plot(map,pathway,color.group = "Type",axis1=1,axis2 = 2,
             plot.title = "Metabolic Pathways-CEED and ACUTE",figure.colors= cols.type,
             show.legend = TRUE,legend.names = names.type,
             pdf.name = "Figures/Figure1/pathway_PCO1&2_all.pdf")

p8<-pco.plot(map,pathway,color.group = "Type",axis1=3,axis2 = 4,
             plot.title = "Metabolic Pathways-CEED and ACUTE",figure.colors= cols.type,
             show.legend = TRUE,legend.names = names.type,
             pdf.name = "Figures/Figure1/pathway_PCO3&4_all.pdf")

#CEED
p9<-pco.plot(map.unc,pathway.unc,color.group = "Type",axis1=1,axis2 = 2,
             plot.title = "Metabolic Pathways-CEED",figure.colors= cols.type,
             show.legend = TRUE,legend.names = names.type,
             pdf.name = "Figures/Figure1/pathway_PCO1&2_CEED.pdf")

p10<-pco.plot(map.unc,pathway.unc,color.group = "Type",axis1=3,axis2 = 4,
              plot.title = "Metabolic Pathways-CEED",figure.colors= cols.type,
              show.legend = TRUE,legend.names = names.type,
              pdf.name = "Figures/Figure1/pathway_PCO3&4_CEED.pdf")
#ACUTE
p11<-pco.plot(map.denver,pathway.denver,color.group = "Type",
              plot.title = "Metabolic Pathways-ACUTE",figure.colors= cols.type,
              show.legend = TRUE,legend.names = names.type,
              pdf.name = "Figures/Figure1/pathway_PCO1&2_ACUTE.pdf")

p12<-pco.plot(map.denver,pathway.denver,color.group = "Type", axis1=3,axis2 = 4,
              plot.title = "Metabolic Pathways-ACUTE",figure.colors= cols.type,
              show.legend = TRUE,legend.names = names.type,
              pdf.name = "Figures/Figure1/pathway_PCO3&4_ACUTE.pdf")

##########Boxplot for phenotypes##########
#BMI
p.bmi<-phenotype.boxplot(map,"BMI","BMI (kg/m2)","CEED and ACUTE")
p.bmi.UNC<-phenotype.boxplot(map.unc,"BMI","BMI (kg/m2)","CEED")
p.bmi.Denver<-phenotype.boxplot(map.denver,"BMI","BMI (kg/m2)","ACUTE")
pdf("Figures/Figure1/bmi.pdf",height = 4, width = 4)
print(p.bmi)
print(p.bmi.UNC)
print(p.bmi.Denver)
dev.off()


# Extended Figure 1 -----------------------------------------------------------------

##########Comparing BMI between UNC and Denver##########
map$Location_name <- ifelse(map$Location == 'Denver','ACUTE','CEED')
map$Type_name <- factor(ifelse(map$Type == 'HC','non-ED',ifelse(map$Type == 'T1','T1','T2')))

compare_locations <- map %>%
  group_by(Type_name) %>%
  rstatix::t_test(BMI ~ Location)  %>%
  adjust_pvalue(method = 'BH') %>%
  add_xy_position() %>% mutate(p.adj.sig = c('ns','***','***'))


plot <- ggplot(map,aes(x=Location_name, y = BMI))+
  geom_boxplot(aes(color=Type_name),outlier.shape=NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3),aes(color=Type_name),alpha=0.5,size=0.5)+
  scale_color_manual(values=cols.type)+
  theme_classic(10)+
  labs(y = "BMI (kg/m2)")+
  facet_wrap(.~Type_name)+
  theme(legend.position = 'none',axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="")+stat_pvalue_manual(compare_locations,label = 'p.adj.sig')

##########UNC phenotypes##########
plot.list<-list()
index<-1
plot.list[[index]]<-plot
index <- index + 1

unc.meta.variables1<-unc.meta.variables[-1]
unc.meta.variables.names1<-unc.meta.variables.names[-1]

for (i in 1:length(unc.meta.variables1)){
  plot<-phenotype.boxplot(map.unc,unc.meta.variables1[i],unc.meta.variables.names1[i],"CEED")
  plot.list[[index]]<-plot
  index<-index+1
}

pdf("Figures/Figure1/phenotype_boxplots.pdf",height = 10,width = 10)
plot_grid(plotlist=plot.list,ncol = 3,nrow=3)
dev.off()


#alpha diversity
shannon<-phenotype.boxplot(map,"shannon_taxonomy","Shannon Diversity Index","Taxonomy_CEED and ACUTE")
shannon.unc<-phenotype.boxplot(map.unc,"shannon_taxonomy","Shannon Diversity Index","Taxonomy_CEED")
shannon.denver<-phenotype.boxplot(map.denver,"shannon_taxonomy","Shannon Diversity Index","Taxonomy_ACUTE")
shannon.pathways<-phenotype.boxplot(map,"shannon_pathways","Shannon Diversity Index","Pathways_CEED and ACUTE")
shannon.unc.pathways<-phenotype.boxplot(map.unc,"shannon_pathways","Shannon Diversity Index","Pathways_CEED")
shannon.denver.pathways<-phenotype.boxplot(map.denver,"shannon_pathways","Shannon Diversity Index","Pathways_ACUTE")

pdf("Figures/Figure1/alpha_diversity.pdf",height = 7,width = 10)
plot_grid(shannon,shannon.unc,shannon.denver,shannon.pathways,shannon.unc.pathways,shannon.denver.pathways,ncol = 3,nrow=2)
dev.off()

pdf("Figures/Figure1/Extended_Figure1.pdf",height = 15,width = 15)
plot_grid(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],
          plot.list[[5]],plot.list[[6]],plot.list[[7]],plot.list[[8]],
          shannon,shannon.unc,shannon.denver,NULL,shannon.pathways,
          shannon.unc.pathways,shannon.denver.pathways,ncol = 4,nrow=4)
dev.off()

# Table 1 -----------------------------------------------------------------
#Missing BMIs:
map$SampleID[is.na(map$BMI)]

tab.summary.cohort <- map %>% group_by(Type,Location) %>%
  summarise(Number = n(),BMI_mean = mean(BMI,na.rm = TRUE),BMI_sd = sd(BMI,na.rm = TRUE),
          Weight_mean = mean(Weight.kg.,na.rm = TRUE),Weight_sd = sd(Weight.kg.,na.rm = TRUE),
          missing_BMI = sum(is.na(BMI)),
          missing_Weight = sum(is.na(Weight.kg.)),
          recovery_days_mean = mean(TimeBetweenRecovery, na.rm = TRUE),recovery_days_sd = sd(TimeBetweenRecovery, na.rm = TRUE),
          BMI_change_mean = mean(BMIchange,na.rm = TRUE), BMI_change_sd = sd(BMIchange,na.rm = TRUE),
          BMI_change_per_day_mean = mean(BMIchangePerDay,na.rm = TRUE), BMI_change_per_day_sd = sd(BMIchangePerDay,na.rm = TRUE),
          weight_change_mean = mean(changeInWeight,na.rm = TRUE), weight_change_sd = sd(changeInWeight,na.rm = TRUE),
          weight_change_per_day_mean = mean(changeInWeightPerDay,na.rm = TRUE), weight_change_per_day_sd = sd(changeInWeightPerDay,na.rm = TRUE),
          TotalFatPercent_mean = mean(TotalFatPercent,na.rm = TRUE), TotalFatPercent_sd = sd(TotalFatPercent, na.rm = TRUE),
          StandardBIA_FatPercent_mean = mean(StandardBIA_FatPercent,na.rm = TRUE), StandardBIA_FatPercent_sd = sd(StandardBIA_FatPercent, na.rm = TRUE),
          Trunc_FatPercent_mean = mean(Trunk..Percent.Fat,na.rm = TRUE), Trunc_FatPercent_sd = sd(Trunk..Percent.Fat, na.rm = TRUE),
          Head_FatPercent_mean = mean(Head..Percent.Fat,na.rm = TRUE), Head_FatPercent_sd = sd(Head..Percent.Fat, na.rm = TRUE),
          Trunc_Fatg_mean = mean(Trunk..Fat.mass..g.,na.rm = TRUE), Trunc_Fatg_sd = sd(Trunk..Fat.mass..g., na.rm = TRUE),
          Head_Fatg_mean = mean(Head..Fat.mass..g.,na.rm = TRUE), Head_Fatg_sd = sd(Head..Fat.mass..g., na.rm = TRUE))


write.table(tab.summary.cohort,"Figures/Figure1/table_summary.txt",sep = '\t',row.names = FALSE)

#Looking at T1 T2 HC overall
tab.summary <- map %>% group_by(Type) %>%
  summarise(BMI_mean = mean(BMI,na.rm = TRUE),BMI_sd = sd(BMI,na.rm = TRUE),
            Weight_mean = mean(Weight.kg.,na.rm = TRUE), Weight_sd = sd(Weight.kg.,na.rm = TRUE))

weight_change <- map %>% filter(Type == 'T1') %>% summarise(mean(changeInWeight,na.rm = TRUE), sd(changeInWeight,na.rm = TRUE))
low_bmi <- map %>% filter(Type == 'T2' & BMI <18.5) %>% group_by(Location) %>% summarise(n())
low_bmi_range <- map %>% filter(Type == 'T2' & BMI <18.5) %>% summarise(min(BMI),max(BMI))

high_bmi <- map %>% filter(Type == 'T2' & BMI >= 18.5) %>% group_by(Location) %>% summarise(n())
high_bmi_range <- map %>% filter(Type == 'T2' & BMI >= 18.5) %>% summarise(min(BMI),max(BMI))


# Supplemental Figure 1 --------------------------------------------------------
p4<-pco.plot(map.unc,species.unc,color.group = "Type", axis1=3,axis2 = 4,
             plot.title = "Taxonomy-CEED",figure.colors = cols.type,
             show.legend = TRUE,legend.names = names.type,
             pdf.name = "Figures/Figure1/Supplemental_Figure1.pdf")

