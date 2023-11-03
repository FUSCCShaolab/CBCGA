rm(list = ls())

library(tidyverse)
library(reshape2)
library(pheatmap)
library(patchwork)
library(ggpubr)
library(plyr)
# 
load("CBCGA.Extended_MergedData_V2.6_220802.Rdata")

info_subtype <- CBCGA_Cohort.Info[, c(1, 10:11, 44)]
info_subtype$`Intrinsic subtype (AIMS)`[info_subtype$`Intrinsic subtype (AIMS)` == "NA"] <- NA
info_subtype <- info_subtype[!is.na(info_subtype$`Intrinsic subtype (PAM50)`), ]
names(info_subtype) <- c("PatientCode", "PAM50 subtypes", "AIMS subtypes", "IHC subtypes")

info_subtype$`PAM50 subtypes` <- factor(info_subtype$`PAM50 subtypes`,
                                        levels = names(Color_PAM50))
info_subtype$`AIMS subtypes` <- factor(info_subtype$`AIMS subtypes`,
                                       levels = names(Color_PAM50))
info_subtype$`IHC subtypes` <- factor(info_subtype$`IHC subtypes`,
                                      levels = names(Color_ClinSub))

# 
fun_confusion <- function(var1, var2){
  confusion1 <- as.data.frame(table(var1, var2))
  names(confusion1) <- c("Var1", "Var2", "Freq")
  # confusion1 <- dcast(data = confusion1,
  #                     Var1~Var2)
  # confusion1 <- column_to_rownames(confusion1, var = "Var1")
  return(confusion1)
}
confu1 <- fun_confusion(info_subtype$`PAM50 subtypes`, info_subtype$`IHC subtypes`)

var1 <- "PAM50 subtypes"
var2 <- "IHC subtypes"
p_1 <- ggplot(confu1)+
  aes(x = Var1, y = Var2)+
  geom_tile(aes(fill = Freq), color = "black")+
  scale_fill_gradient(low = "white", high = "darkred")+
  scale_y_discrete(limits = rev(names(Color_ClinSub)))+#
  geom_text(aes(label = Freq))+
  labs(x = var1, y = var2)+
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks.length = unit(0, "cm"),
        legend.position = "none");p_1

p_var1 <- ggplot(confu1)+
  aes(x = Var1, y = Freq, fill = Var2)+
  geom_bar(width = 1, color = "black",
           stat = "identity", position = "fill")+
  scale_fill_manual(values = Color_ClinSub)+
  scale_y_continuous(expand = c(0, 0), n.breaks = 6)+
  labs(fill = "")+
  theme_classic()+
  theme(
    legend.position = "left",
    axis.line.x = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.length.x = unit(0, "cm"));p_var1
p_var2 <- ggplot(confu1)+
  aes(x = Var2, y = Freq)+
  geom_bar(aes(fill = Var1),width = 1, color = "black",
           stat = "identity", position = "fill")+
  scale_fill_manual(values = Color_PAM50)+
  scale_y_continuous(expand = c(0, 0), n.breaks = 6)+
  scale_x_discrete(limits = rev(names(Color_ClinSub)))+#
  coord_flip()+
  labs(fill = "")+
  theme_classic()+
  theme(
    axis.line.y = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "cm"));p_var2

ggarrange(p_var1, NULL,
          p_1, p_var2,
          nrow = 2, ncol = 2,
          heights = c(.5,1), 
          widths = c(1,.5))
ggsave("Results/CBCGA_confusion_table_PAM50-IHC.pdf", width = 8, height = 6)
openxlsx::write.xlsx(confu1,file = "ED_Fig1B.xlsx")

# --------
confu1 <- fun_confusion(info_subtype$`PAM50 subtypes`, info_subtype$`AIMS subtypes`)
var1 <- "PAM50 subtypes"
var2 <- "AIMS subtypes"
p_1 <- ggplot(confu1)+
  aes(x = Var1, y = Var2)+
  geom_tile(aes(fill = Freq), color = "black")+
  scale_fill_gradient(low = "white", high = "darkred")+
  scale_y_discrete(limits = rev(names(Color_PAM50)))+#
  geom_text(aes(label = Freq))+
  labs(x = var1, y = var2)+
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks.length = unit(0, "cm"),
        legend.position = "none");p_1

p_var1 <- ggplot(confu1)+
  aes(x = Var1, y = Freq)+
  geom_bar(aes(fill = Var2),width = 1, color = "black",
           stat = "identity", position = "fill")+
  scale_fill_manual(values = Color_PAM50)+
  scale_y_continuous(expand = c(0, 0), n.breaks = 6)+
  labs(fill = "")+
  theme_classic()+
  theme(
    legend.position = "left",
    axis.line.x = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.length.x = unit(0, "cm"));p_var1
p_var2 <- ggplot(confu1)+
  aes(x = Var2, y = Freq)+
  geom_bar(aes(fill = Var1),width = 1, color = "black",
           stat = "identity", position = "fill")+
  scale_fill_manual(values = Color_PAM50)+
  scale_y_continuous(expand = c(0, 0), n.breaks = 6)+
  scale_x_discrete(limits = rev(names(Color_PAM50)))+#
  coord_flip()+
  labs(fill = "")+
  theme_classic()+
  theme(
    # legend.position = "none",
    axis.line.y = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "cm"));p_var2

ggarrange(p_var1, NULL,
          p_1, p_var2,
          nrow = 2, ncol = 2,
          heights = c(.5,1), 
          widths = c(1,.5))
ggsave("Results/CBCGA_confusion_table_PAM50-AIMS.pdf", width = 8, height = 6)
openxlsx::write.xlsx(confu1,file = "ED_Fig1C.xlsx")
