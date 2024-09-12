
library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)

setwd("T:\\TGB\\LSS\\TGU\\Users\\Tomas\\2024_VolcanoPlot_Mozhgan")

table <- read_excel("Volcano plot.xlsx",
                    "BP QUANTIFIABLE with ID",
                    skip=1,
)

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

ggplot(table, aes(`Log2(Sen/Pro)`, -log(`P-Value Pro vs Sen`), col=SenVPro_EffectClass)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values = c("grey", "blue", "red"), na.translate = FALSE) +
  geom_text_repel(aes(label = SenVPro_label), show.legend=FALSE) +
  xlab("Log2 Fold Change") +
  ylab("Negative Log p-Value") +
  labs(title = "Senescent vs Proliferating", col="Effect Direction")

ggplot(table, aes(`Log2(Sen/Qui)`, -log(`P-Value  Qui vs Sen`), col=SenVQui_EffectClass)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("grey", "blue", "red"), na.translate = FALSE) +
  geom_text_repel(aes(label = SenVQui_label), show.legend=FALSE) + 
  xlab("Log2 Fold Change") + 
  ylab("Negative Log p-Value") +
  labs(title = "Senescent vs Quiescent", col="Effect Direction")

ggplot(table, aes(`Log2(Pro/Qui)`, -log(`P-Value  Pro vs Qui`), col=ProVQui_EffectClass)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("grey", "blue", "red"), na.translate = FALSE) +
  geom_text_repel(aes(label = ProVQui_label), show.legend=FALSE) + 
  xlab("Log2 Fold Change") +
  ylab("Negative Log p-Value") +
  labs(title = "Proliferating vs Quiescent", col="Effect Direction")
