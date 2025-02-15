
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(RColorBrewer)
library(pheatmap)

setwd("T:\\TGB\\LSS\\TGU\\Users\\Tomas\\2024_TopDown_Mozhgan")

table <- read_excel("HeatMap_normalized.xlsx", 
                    "QUANTIFIABLE with ID")

table2 <- pivot_longer(table, cols=c(`Pro 1`, 
                                     `Pro 2`, 
                                     `Pro 3`, 
                                     `Pro 4`, 
                                     `Qui 1`, 
                                     `Qui 2`, 
                                     `Qui 3`, 
                                     `Qui 4`, 
                                     `SEN 1`, 
                                     `SEN 2`, 
                                     `SEN 3`, 
                                     `SEN 4`),
                       names_to="cond.rep",
                       values_to="intensity")

table2[c('Cond', 'Rep')] <- str_split_fixed(table2$cond.rep, ' ', 2)

table2$intensity[is.na(table2$intensity)] = 0.0

table2$Cond <- str_to_title(table2$Cond)

mtrx <- with(table2, tapply(intensity, list(Cond, mass), mean))

mtrx2 <- t(mtrx)

mtrx3 <- with(table2, tapply(intensity, list(cond.rep, mass), mean))

mtrx4 <- t(mtrx3)

mtrx4.dist <- dist(mtrx4)

mtrx4.clusters <- hclust(mtrx4.dist, method="ward.D2")

clusterids <- cutree(mtrx4.clusters, k=64)

df <- as.data.frame(mtrx4)

df.withOrd <- cbind(df, ordering=mtrx4.clusters$order, clusterid=clusterids)
df.ordered <- sort_by(df.withOrd, y=~ordering)

res <- pheatmap(mtrx4,
         kmeans_k=NA,
         cluster_rows=mtrx4.clusters,
         cluster_cols=FALSE,
         show_rownames=FALSE,
         scale="row",
         color=colorRampPalette(c("blue", "white", "red"))(10000),
         cellheight=2
)

write.csv(df.ordered, "heatmaptable_inorder.csv")



