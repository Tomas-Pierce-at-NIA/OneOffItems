
library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(stringr)
library(tidyr)
library(RColorBrewer)
library(pheatmap)

setwd("T:\\TGB\\LSS\\TGU\\Users\\Tomas\\2024_TopDown_Mozhgan")

table <- read_excel("HeatMap.xlsx", 
                    "QUANTIFIABLE with ID",
                    skip=1)

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

pheatmap(mtrx4,
         kmeans_k=NA,
         cluster_rows=TRUE,
         cluster_cols=FALSE,
         show_rownames=FALSE,
         scale="row",
         color=colorRampPalette(c("blue", "white", "red"))(10000),
         cellheight=2
)

res <- pheatmap(mtrx4,
                kmeans_k=NA,
                cluster_rows=TRUE,
                cluster_cols=FALSE,
                show_rownames=FALSE,
                scale="row",
                color=colorRampPalette(c("blue", "white", "red"))(10000),
                cellheight=2
)

mtrx4.clust <- cbind(mtrx4,
                     clust = cutree(res$tree_row,
                                    h=4))

write.table(mtrx4.clust, "pheatclust.txt", row.names=TRUE, col.names=TRUE)
