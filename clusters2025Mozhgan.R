library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("T:\\TGB\\LSS\\TGU\\Users\\Tomas\\2024_TopDown_Mozhgan")

table1 <- read_excel("pathway_protein_lists.xlsx",
                    sheet="pathway_protein_lists (002)")

table2 <- read_excel("pathway_protein_lists.xlsx",
                     sheet="total identified proteins")

senvqui <- table1$SenVQui
senvqui <- na.omit(senvqui)
senvqui.enrz <- bitr(senvqui,
                     "SYMBOL",
                     "ENTREZID",
                     "org.Hs.eg.db")
senvpro <- na.omit(table1$SenVPro)
senvpro.enrz <- bitr(senvpro,
                     "SYMBOL",
                     "ENTREZID",
                     "org.Hs.eg.db")

provqui <- na.omit(table1$ProVQui)
provqui.enrz <- bitr(provqui,
                     "SYMBOL",
                     "ENTREZID",
                     "org.Hs.eg.db")

senexcl <- na.omit(table1$Sen_Exclusive)
senexcl.enrz <- bitr(provqui,
                     "SYMBOL",
                     "ENTREZID",
                     "org.Hs.eg.db")

verse <- table2$gene
verse.enrz <- bitr(verse,
                   "SYMBOL",
                   "ENTREZID",
                   "org.Hs.eg.db")

# done correctly

cl_provqui1 <- enrichGO(provqui.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="MF",
                        universe=verse.enrz$ENTREZID)

dotplot(cl_provqui1,
        title="Proliferating vs Quiescent Molecular Function")

cl_provqui2 <- enrichGO(provqui.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="CC",
                        universe=verse.enrz$ENTREZID)

dotplot(cl_provqui2,
        title="Proliferating vs Quiescent Cellular Compartment")

cl_provqui3 <- enrichGO(provqui.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="BP",
                        universe=verse.enrz$ENTREZID)

dotplot(cl_provqui3,
        title="Proliferating vs Quiescent Biological Process")

cl_senexcl1 <- enrichGO(senexcl.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="MF",
                        universe=verse.enrz$ENTREZID)

dotplot(cl_senexcl1,
        title="Senescent Exclusive Molecular Function")

cl_senexcl2 <- enrichGO(senexcl.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="CC",
                        universe=verse.enrz$ENTREZID)

dotplot(cl_senexcl2,
        title="Senescent Exclusive Cellular Compartment")

cl_senexcl3 <- enrichGO(senexcl.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="BP",
                        universe=verse.enrz$ENTREZID)

dotplot(cl_senexcl3,
        title="Senescent Exclusive Biological Process")

cl_senvpro1 <- enrichGO(senvpro.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="MF",
                        universe=verse.enrz$ENTREZID)

dotplot(cl_senvpro1,
        title="Senescent versus Proliferating Molecular Function")

cl_senvpro2 <- enrichGO(senvpro.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="CC",
                        universe=verse.enrz$ENTREZID)

dotplot(cl_senvpro2,
        title="Senescent versus Proliferating Cellular Compartment")

cl_senvpro3 <- enrichGO(senvpro.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="BP",
                        universe=verse.enrz$ENTREZID)

dotplot(cl_senvpro3,
        title="Senescent versus Profilerating Biological Process")

cl_senvqui1 <- enrichGO(senvqui.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="MF",
                        universe=verse.enrz$ENTREZID)

dotplot(cl_senvqui1,
        title="Senescent versus Quiescent Molecular Function")

cl_senvqui2 <- enrichGO(senvqui.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="CC",
                        universe=verse.enrz$ENTREZID)

dotplot(cl_senvqui2,
        title="Senescent versus Quiescent Cellular Compartment")

cl_senvqui3 <- enrichGO(senvqui.enrz$ENTREZID,
                         "org.Hs.eg.db",
                         ont="BP",
                         universe=verse.enrz$ENTREZID)

dotplot(cl_senvqui3,
        title="Senescent versus Quiescent Biological Process")


# not corrected


cl_provqui1 <- enrichGO(provqui.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="MF",
                        universe=verse.enrz$ENTREZID,
                        pAdjustMethod="none")

dotplot(cl_provqui1,
        title="Proliferating vs Quiescent Molecular Function uncorrected")

cl_provqui2 <- enrichGO(provqui.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="CC",
                        universe=verse.enrz$ENTREZID,
                        pAdjustMethod="none")

dotplot(cl_provqui2,
        title="Proliferating vs Quiescent Cellular Compartment uncorrected")

cl_provqui3 <- enrichGO(provqui.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="BP",
                        universe=verse.enrz$ENTREZID,
                        pAdjustMethod="none")

dotplot(cl_provqui3,
        title="Proliferating vs Quiescent Biological Process uncorrected")

cl_senexcl1 <- enrichGO(senexcl.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="MF",
                        universe=verse.enrz$ENTREZID,
                        pAdjustMethod="none")

dotplot(cl_senexcl1,
        title="Senescent Exclusive Molecular Function uncorrected")

cl_senexcl2 <- enrichGO(senexcl.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="CC",
                        universe=verse.enrz$ENTREZID,
                        pAdjustMethod="none")

dotplot(cl_senexcl2,
        title="Senescent Exclusive Cellular Compartment uncorrected")

cl_senexcl3 <- enrichGO(senexcl.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="BP",
                        universe=verse.enrz$ENTREZID,
                        pAdjustMethod="none")

dotplot(cl_senexcl3,
        title="Senescent Exclusive Biological Process uncorrected")

cl_senvpro1 <- enrichGO(senvpro.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="MF",
                        universe=verse.enrz$ENTREZID,
                        pAdjustMethod="none")

dotplot(cl_senvpro1,
        title="Senescent versus Proliferating Molecular Function uncorrected")

cl_senvpro2 <- enrichGO(senvpro.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="CC",
                        universe=verse.enrz$ENTREZID,
                        pAdjustMethod="none")

dotplot(cl_senvpro2,
        title="Senescent versus Proliferating Cellular Compartment uncorrected")

cl_senvpro3 <- enrichGO(senvpro.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="BP",
                        universe=verse.enrz$ENTREZID,
                        pAdjustMethod="none")

dotplot(cl_senvpro3,
        title="Senescent versus Profilerating Biological Process uncorrected")

cl_senvqui1 <- enrichGO(senvqui.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="MF",
                        universe=verse.enrz$ENTREZID,
                        pAdjustMethod="none")

dotplot(cl_senvqui1,
        title="Senescent versus Quiescent Molecular Function uncorrected")

cl_senvqui2 <- enrichGO(senvqui.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="CC",
                        universe=verse.enrz$ENTREZID,
                        pAdjustMethod="none")

dotplot(cl_senvqui2,
        title="Senescent versus Quiescent Cellular Compartment uncorrected")

cl_senvqui3 <- enrichGO(senvqui.enrz$ENTREZID,
                        "org.Hs.eg.db",
                        ont="BP",
                        universe=verse.enrz$ENTREZID,
                        pAdjustMethod="none")

dotplot(cl_senvqui3,
        title="Senescent versus Quiescent Biological Process uncorrected")
