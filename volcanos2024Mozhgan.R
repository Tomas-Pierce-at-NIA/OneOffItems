
library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(stringr)

setwd("T:\\TGB\\LSS\\TGU\\Users\\Tomas\\2024_TopDown_Mozhgan")

table <- read_excel("20250121-Volcano plots.xlsx",
                    "BP QUANTIFIABLE with ID",
                    skip=1,
)

#table[c('Gene.Name', 'PTM')] <- str_split_fixed(table$GENE, ' ', 2)

colnames(table)[3] <- "PTM"

table$PTM[is.na(table$PTM)] = "NOTLISTED"

table <- table %>%
  mutate(SenVPro_EffectClass = ifelse(`P-Value Pro vs Sen` < 0.05,
                                      ifelse(`Log2(Sen/Pro)` > 0.0,
                                             "Senescence enhanced",
                                             "Senescence depleted"),
                                      "no effect"))

table <- table %>%
  mutate(SenVQui_EffectClass = ifelse(`P-Value  Qui vs Sen` < 0.05,
                                      ifelse(`Log2(Sen/Qui)` > 0.0,
                                             "Senescence enhanced",
                                             "Senescence depleted"),
                                      "no effect"))

table <- table %>%
  mutate(ProVQui_EffectClass = ifelse(`P-Value  Pro vs Qui` < 0.05,
                                      ifelse(`Log2(Pro/Qui)` > 0.0,
                                             "Proliferation enhanced",
                                             "Proliferation depleted"),
                                      "no effect"))

table <- table %>%
  mutate(SenVPro_label = ifelse(SenVPro_EffectClass == "no effect",
                                "",
                                GENE))

table <- table %>%
  mutate(SenVQui_label = ifelse(SenVQui_EffectClass == "no effect",
                                "",
                                GENE))

table <- table %>%
  mutate(ProVQui_label = ifelse(ProVQui_EffectClass == "no effect",
                                "",
                                GENE))

#table$PTM[table$PTM == "(P1)"] <- "Proteoform 1"
#table$PTM[table$PTM == "(P2)"] <- "Proteoform 2"
#table$PTM[table$PTM == "(P3)"] <- "Proteoform 3"


pdf("SenVPro.pdf", width=36, height=36)

ggplot(table, aes(`Log2(Sen/Pro)`, -log(`P-Value Pro vs Sen`), col=SenVPro_EffectClass,
                  shape=PTM)) +
  geom_point(size=8.0) + 
  theme_minimal() +
  scale_shape_manual(
    breaks = c("Proteoform 1", 
               "Proteoform 2", 
               "Proteoform 3"),
    values = c("NOTLISTED" = 16, 
               "Proteoform 1" = 15, 
               "Proteoform 2" = 17, 
               "Proteoform 3" = 3)
  ) + 
  scale_color_manual(values = c("grey", "blue", "red"), na.translate = FALSE) +
  geom_text_repel(aes(label = SenVPro_label), 
                  show.legend=FALSE,
                  size=18,
                  point.padding=unit(0.5, "inches"),
                  box.padding=unit(0.25, "inches"),
                  force=2,
                  max.time=10,
                  max.iter=100000) +
  xlab("Log2 Fold Change") +
  ylab("Negative Log p-Value") +
  labs(title = "Senescent vs Proliferating", col="Effect Direction") +
  geom_hline(yintercept = -log(0.05)) + 
  theme(axis.title = element_text(size=55),
        axis.text = element_text(size=45),
        plot.title = element_text(size=60),
        legend.title = element_text(size=36),
        legend.text = element_text(size=36)) +
  xlim(-6, 6) +
  ylim(0, 12.5)

dev.off()

pdf("SenVPro_unlabeled.pdf", width=36, height=36)

ggplot(table, aes(`Log2(Sen/Pro)`, -log(`P-Value Pro vs Sen`), col=SenVPro_EffectClass,
                  shape=PTM)) +
  geom_point(size=8.0) + 
  theme_minimal() +
  scale_shape_manual(
    breaks = c("Proteoform 1", 
               "Proteoform 2", 
               "Proteoform 3"),
    values = c("NOTLISTED" = 16, 
               "Proteoform 1" = 15, 
               "Proteoform 2" = 17, 
               "Proteoform 3" = 3)
  ) + 
  scale_color_manual(values = c("grey", "blue", "red"), na.translate = FALSE) +
  # geom_text_repel(aes(label = SenVPro_label), 
  #                 show.legend=FALSE,
  #                 size=18,
  #                 point.padding=unit(0.5, "inches"),
  #                 box.padding=unit(0.25, "inches"),
  #                 force=2,
  #                 max.time=10,
  #                 max.iter=100000) +
  xlab("Log2 Fold Change") +
  ylab("Negative Log p-Value") +
  labs(title = "Senescent vs Proliferating", col="Effect Direction") +
  geom_hline(yintercept = -log(0.05)) + 
  theme(axis.title = element_text(size=55),
        axis.text = element_text(size=45),
        plot.title = element_text(size=60),
        legend.title = element_text(size=36),
        legend.text = element_text(size=36)) +
  xlim(-6, 6) +
  ylim(0, 12.5)


dev.off()

pdf("QuiVSen.pdf", width=36, height=36)

ggplot(table, aes(`Log2(Sen/Qui)`, -log(`P-Value  Qui vs Sen`), col=SenVQui_EffectClass,
                  shape=PTM)) +
  geom_point(size=8.0) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "blue", "red"), na.translate = FALSE) +
  scale_shape_manual(
    breaks = c("Proteoform 1", 
               "Proteoform 2", 
               "Proteoform 3"),
    values = c("NOTLISTED" = 16, 
               "Proteoform 1" = 15, 
               "Proteoform 2" = 17, 
               "Proteoform 3" = 3)
  ) +
  geom_text_repel(aes(label = SenVQui_label), 
                  show.legend=FALSE,
                  size=18,
                  point.padding=unit(0.5, "inches"),
                  box.padding=unit(0.25, "inches"),
                  force=2,
                  max.time=10,
                  max.iter=100000) + 
  xlab("Log2 Fold Change") + 
  ylab("Negative Log p-Value") +
  labs(title = "Senescent vs Quiescent", col="Effect Direction") +
  geom_hline(yintercept = -log(0.05)) + 
  theme(axis.title = element_text(size=50),
        axis.text = element_text(size=40),
        plot.title = element_text(size=60),
        legend.title = element_text(size=36),
        legend.text = element_text(size=36)) +
  xlim(-6, 6) +
  ylim(0, 12.5)

dev.off()

pdf("QuiVSen_unlabeled.pdf", width=36, height=36)

ggplot(table, aes(`Log2(Sen/Qui)`, -log(`P-Value  Qui vs Sen`), col=SenVQui_EffectClass,
                  shape=PTM)) +
  geom_point(size=8.0) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "blue", "red"), na.translate = FALSE) +
  scale_shape_manual(
    breaks = c("Proteoform 1", 
               "Proteoform 2", 
               "Proteoform 3"),
    values = c("NOTLISTED" = 16, 
               "Proteoform 1" = 15, 
               "Proteoform 2" = 17, 
               "Proteoform 3" = 3)
  ) +
  # geom_text_repel(aes(label = SenVQui_label), 
  #                 show.legend=FALSE,
  #                 size=14,
  #                 point.padding=unit(0.5, "inches"),
  #                 box.padding=unit(0.25, "inches"),
  #                 force=2,
  #                 max.time=10,
  #                 max.iter=100000) + 
  xlab("Log2 Fold Change") + 
  ylab("Negative Log p-Value") +
  labs(title = "Senescent vs Quiescent", col="Effect Direction") +
  geom_hline(yintercept = -log(0.05)) + 
  theme(axis.title = element_text(size=50),
        axis.text = element_text(size=40),
        plot.title = element_text(size=60),
        legend.title = element_text(size=36),
        legend.text = element_text(size=36)) +
  xlim(-6, 6) +
  ylim(0, 12.5)

dev.off()


pdf("ProVQui.pdf", width=36, height=36)

ggplot(table, aes(`Log2(Pro/Qui)`, -log(`P-Value  Pro vs Qui`), col=ProVQui_EffectClass,
                  shape=PTM)) +
  geom_point(size=8.0) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "blue", "red"), na.translate = FALSE) +
  geom_text_repel(aes(label = ProVQui_label), 
                  show.legend=FALSE,
                  size=18,
                  point.padding=unit(0.5, "inches"),
                  box.padding=unit(0.25, "inches"),
                  force=2,
                  max.time=10,
                  max.iter=100000) + 
  scale_shape_manual(
    breaks = c("Proteoform 1", 
               "Proteoform 2", 
               "Proteoform 3"),
    values = c("NOTLISTED" = 16, 
               "Proteoform 1" = 15, 
               "Proteoform 2" = 17, 
               "Proteoform 3" = 3)
  ) +
  xlab("Log2 Fold Change") +
  ylab("Negative Log p-Value") +
  labs(title = "Proliferating vs Quiescent", col="Effect Direction") + 
  geom_hline(yintercept = -log(0.05)) + 
  theme(axis.title = element_text(size=55),
        axis.text = element_text(size=45),
        plot.title = element_text(size=60),
        legend.title = element_text(size=36),
        legend.text = element_text(size=36)) +
  xlim(-6, 6) + 
  ylim(0, 12.5)

dev.off()



pdf("ProVQui_unlabeled.pdf", width=36, height=36)

ggplot(table, aes(`Log2(Pro/Qui)`, -log(`P-Value  Pro vs Qui`), col=ProVQui_EffectClass,
                  shape=PTM)) +
  geom_point(size=8.0) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "blue", "red"), na.translate = FALSE) +
  # geom_text_repel(aes(label = ProVQui_label), 
  #                 show.legend=FALSE,
  #                 size=14,
  #                 point.padding=unit(0.5, "inches"),
  #                 box.padding=unit(0.25, "inches"),
  #                 force=2,
  #                 max.time=10,
  #                 max.iter=100000) + 
  scale_shape_manual(
    breaks = c("Proteoform 1", 
               "Proteoform 2", 
               "Proteoform 3"),
    values = c("NOTLISTED" = 16, 
               "Proteoform 1" = 15, 
               "Proteoform 2" = 17, 
               "Proteoform 3" = 3)
  ) +
  xlab("Log2 Fold Change") +
  ylab("Negative Log p-Value") +
  labs(title = "Proliferating vs Quiescent", col="Effect Direction") + 
  geom_hline(yintercept = -log(0.05)) + 
  theme(axis.title = element_text(size=55),
        axis.text = element_text(size=45),
        plot.title = element_text(size=60),
        legend.title = element_text(size=36),
        legend.text = element_text(size=36)) +
  xlim(-6, 6) + 
  ylim(0, 12.5)

dev.off()

