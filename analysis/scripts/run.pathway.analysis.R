
moduleDir = dirname(getwd())
pipeRoot = dirname(moduleDir)

pathway.path = file.path(pipeRoot, "input","pathway", "humanN2_pathabundance_cpm.tsv")
meta.path = file.path(pipeRoot, "meta_master.txt")

outputDir=file.path(moduleDir, "output")

source("../resources/analysis.R")

run.analysis.all("pathway",
                 pathway.path,
                 meta.path,
                 outputDir,
                 is.taxonomy = FALSE)
