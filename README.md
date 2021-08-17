# Anorexia_Analyses2021
This repository includes analyses for the paper **Gut Microbes are Associated with Weight Trajectory in Anorexia Nervosa Therapy**.

## Steps to run the analyses:
1. For each taxonomic level as well as metabolic pathways, we compare the abundances of microbial features between the non-ED, T1, and T2 groups. We also study the associations between the gut microbiome and BMI, weight, BMI change per dat, and other variables. 

Run the following codes in order when you set your work directory at the Anorexia_Analyses2021:

```{r,eval=F, echo=F}
R/run.script.R
R/run.analysis.R
```

2. To generate the figures in the manuscript, run the following codes in order:

```{r,eval=F, echo=F}
R/PrepareDataForFigures.R
R/Figure1.R
R/Figure2.R
R/Figure3.R
R/Figure4.R
```

