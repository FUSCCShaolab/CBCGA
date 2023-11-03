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
col <- c(rep(c("#beced2"),5),rep(c("#6281a5"),19),rep(c("#ca9d9d"),11),rep(c("#b33e38"),2))

Modal_Res_Ord <- Modal_Res[order(Modal_Res$C_index,decreasing = F),]
Modal_Res_Ord$Abbr <- factor(Modal_Res_Ord$Abbr, levels = c(Modal_Res_Ord$Abbr))

ggplot(Modal_Res_Ord,aes(x=Abbr,y=C_index)) +
  geom_bar(stat='identity',fill = col,width = 0.7) +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0.4, position = position_dodge(0.9))+
  theme_base() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+ 
  scale_y_continuous(breaks=seq(0, 1, 0.2))

