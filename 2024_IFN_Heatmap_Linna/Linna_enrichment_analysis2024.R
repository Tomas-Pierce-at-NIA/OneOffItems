
library(enrichR)
library(readxl)
library(clusterProfiler)

setwd("T:/TGB/LSS/TGU/Users/Linna/2022_0602_monocyte_sen_antibodies_list")

ifn_labels <- read_excel("Copy of Proteomic analysis.xlsx",
                         "IFN Heatmap Data")
websiteLive <- getOption("enrichR.live")
setEnrichrSite("Enrichr")

dbl <- listEnrichrDbs()

dbn <- c("GO_Molecular_Function_2017", "GO_Cellular_Component_2017", "GO_Biological_Process_2017")

enriched <- enrichr(ifn_labels$Genes, dbn)

plotEnrich(enriched[[3]], showTerms=20, numChar=40, y = "Count", orderBy="P.value")
plotEnrich(enriched[[2]], showTerms=20, numChar=40, y = "Count", orderBy="P.value")
plotEnrich(enriched[[1]], showTerms=20, numChar=40, y = "Count", orderBy="P.value")
