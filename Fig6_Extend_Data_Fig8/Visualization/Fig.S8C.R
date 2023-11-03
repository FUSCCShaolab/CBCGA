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


Feat_imp <- read.csv("Feature_importance.csv")

ggplot(Feat_imp,aes(x = delta_C, y =reorder(Feat_New,delta_C,decreasing = FALSE), fill = delta_C)) +
  geom_point(size = 5, shape = 21, colour = "grey")  +
  scale_fill_gradient(low = "white",high = "#6281a5") +
  theme_bw() + 
  scale_x_continuous(breaks=seq(0, 0.2, 0.05)) +
  ylab("") +
  xlab("Feature Importance Score")
