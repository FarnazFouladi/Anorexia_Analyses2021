
args = commandArgs(trailingOnly=TRUE)
taxanomicLevel = args[1]

moduleDir = dirname(getwd())
pipeRoot = dirname(moduleDir)
otu.file.name = paste0("anorexia2020Sep17_taxaCount_", taxanomicLevel, ".tsv")
otu.path = file.path(pipeRoot, "input","taxaCounts", otu.file.name)
meta.path = file.path(pipeRoot, "meta_master.txt")

outputDir=file.path(moduleDir, "output")

source("../resources/analysis.R")

run.analysis.all(taxanomicLevel,
                 otu.path,
                 meta.path,
                 outputDir,
                 is.taxonomy = TRUE)
