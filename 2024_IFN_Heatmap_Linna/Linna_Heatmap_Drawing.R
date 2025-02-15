
# September 4 2024
# Tomas Pierce
# Linna Cui



setwd("T:/TGB/LSS/TGU/Users/Linna/2022_0602_monocyte_sen_antibodies_list")

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdendro)

replicates <- read_excel("Copy of Proteomic analysis.xlsx",
                     "Replicate Data")
ifn_labels <- read_excel("Copy of Proteomic analysis.xlsx",
                         "IFN Heatmap Data")
combo <- replicates %>% inner_join(ifn_labels, join_by(PG.Genes == Genes))

focused <- combo %>%
  select(PG.Genes, MCtl1, MCtl2, MCtl3, MCtl4, MCtl5, MSen1, MSen2, MSen3, MSen4, MSen5,
         CategoryType1, CategoryType2, CategoryType3, CategoryType4, CategoryType5, 
         CategoryType6, CategoryType7, CategoryType8, CategoryType9, CategoryType10,
         CategoryType11, CategoryType12, CategoryType13, CategoryType14, CategoryType15,
         CategoryType16, CategoryType17, CategoryType18, CategoryType19, IFNPathwayType)

focused$inIFNPath1 <- agrepl("1", focused$IFNPathwayType, 0.0)
focused$inIFNPath2 <- agrepl("2", focused$IFNPathwayType, 0.0)
focused$inIFNPath3 <- agrepl("3", focused$IFNPathwayType, 0.0)

focused$inCat1 <- !is.na(focused$CategoryType1)

tallfocused <- focused %>% pivot_longer(
  cols = c(MCtl1, MCtl2, MCtl3, MCtl4, MCtl5, MSen1, MSen2, MSen3, MSen4, MSen5),
  names_to="Condition",
  values_to="Value"
) %>% pivot_longer(
  cols = c(CategoryType1, CategoryType2, CategoryType3, CategoryType4, CategoryType5, 
           CategoryType6, CategoryType7, CategoryType8, CategoryType9, CategoryType10,
           CategoryType11, CategoryType12, CategoryType13, CategoryType14, CategoryType15,
           CategoryType16, CategoryType17, CategoryType18, CategoryType19),
  names_to="_CategoryNames",
  values_to="Category"
) %>% na.omit()

tallfocused$isCtrl <- agrepl("Ctl", tallfocused$Condition, 0.0)

tallfocused$WhichCond = ifelse(tallfocused$isCtrl, "Control", "Senescent")

ggplot(tallfocused, aes(WhichCond, PG.Genes, fill=Value)) + 
  geom_tile() + 
  xlab("Condition") +
  ylab("Gene")

ggplot(tallfocused, aes(Condition, PG.Genes, fill=Value)) + 
  geom_tile() + 
  xlab("Replicate") + 
  ylab("Gene")


