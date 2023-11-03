rm(list=ls())

library(stringr)
library(dplyr)
library(glmnet)
library(pROC)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(PMCMRplus)
library(pheatmap)
library(ggthemes)
library(export)
library(survcomp)
library(Hmisc)
library(survminer)

Modal_Res <- read.csv("r0.6_knum5_rs111.csv")
Modal_Res_Ord <- Modal_Res[order(Modal_Res$C_index,decreasing = F),]
Modal_Res_Ord$Abbr <- factor(Modal_Res_Ord$Abbr, levels = c(Modal_Res_Ord$Abbr))

Modal_Res_Ord$Dimension_All <- nchar(as.character(Modal_Res_Ord$Abbr))
Modal_Res_Ord$Dimension_All <- as.numeric(Modal_Res_Ord$Dimension_All)
Modal_Res_Ord$Dimension_Comb <- case_when(Modal_Res_Ord$Dimension_All==1~"Single",
                                          Modal_Res_Ord$Dimension_All==2|Modal_Res_Ord$Dimension_All==3~"Mid",
                                          Modal_Res_Ord$Dimension_All>3 ~"Multi")

Modal_Res_Ord$Dimension_Comb <- factor(Modal_Res_Ord$Dimension_Comb, levels = c("Single", "Mid", "Multi"))
pairwise.wilcox.test(Modal_Res_Ord$C_index,Modal_Res_Ord$Dimension_Comb, p.adjust.method = "fdr",alternative = c("two.sided"))

ggviolin(Modal_Res_Ord, x="Dimension_Comb", y="C_index", fill = "Dimension_Comb",
         palette =  c(colorRampPalette(c("#dfe4ec","#6583a7"))(3)), # "#bec9d9",
         add = "boxplot", add.params = list(fill="white", width = 0.1), width = 1)+
  stat_compare_means(label.y = 1) +
  scale_y_continuous(breaks = c(0.5,0.6,0.7,0.8,0.9,1))

pairwise.wilcox.test(Modal_Res_Ord$C_index,Modal_Res_Ord$Dimension_Comb, p.adjust.method = "fdr")
