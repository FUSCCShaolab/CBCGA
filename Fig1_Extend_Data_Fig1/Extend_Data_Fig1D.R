rm(list = ls())

library(tidyverse)
library(clusterProfiler)
library(DescTools)
library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)
library(xlsx)

# 
load("CBCGA.Extended_MergedData_V2.6_220802.Rdata")
load("Figure1_aimpathway.Rdata")

info_subtype <- CBCGA_Cohort.Info[, c(1, 10:11, 44)]
info_subtype$`Intrinsic subtype (AIMS)`[info_subtype$`Intrinsic subtype (AIMS)` == "NA"] <- NA
names(info_subtype) <- c("PatientCode", "PAM50 subtypes", "AIMS subtypes", "IHC subtypes")

if(TRUE){
  data_protein <- CBCGA.Extended_PRO_normalized
  info_pro <- CBCGA.Extended_PRO_label
  
  id_pro <- names(data_protein)
  id_pro <- id_pro[str_detect(id_pro, "_T$")]
  table(duplicated(str_extract(id_pro, "[A-Z]*_")))
  
  info_subtype$proID <- id_pro[match(info_subtype$PatientCode, str_extract(id_pro, "[A-Z]{4,4}"))]
  info_subtype <- info_subtype[!is.na(info_subtype$proID), ]
  info_subtype <- info_subtype[!is.na(info_subtype$`PAM50 subtypes`), ]
  
  data_protein <- data_protein[info_subtype$proID]
  data_protein <- data_protein[match(info_subtype$proID, names(data_protein))]
  identical(info_subtype$proID, names(data_protein))
  
  info_subtype$PAM50 <- info_subtype$`PAM50 subtypes`
  
  contrast_group <- unique(info_subtype$PAM50)
  
  con <- contrast_group
  
  fun_sig <- function(dat, logFC = 1, FDR = 0.05){
    res <- ifelse(dat$FDR < FDR & abs(dat$logFC)> logFC,
                  ifelse(dat$logFC > 0, "Up", "Down"),
                  "Not")
    return(res)
  }
  
  res_diff <- map(con, function(gg){
    con_info <- info_subtype
    con_info$PAM50 <- ifelse(con_info$PAM50 == gg, gg, "Other")
    con_info$PAM50 <- factor(con_info$PAM50, levels = c(gg, "Other"))
    dat <- data_protein[names(data_protein)%in%con_info$proID]
    
    DEG <- map(seq_along(rownames(dat)), function(x){
      symbol <- rownames(dat)[x]
      dat <- unlist(dat[x, ])
      group <- con_info$PAM50
      test <- wilcox.test(dat~group)
      pvalue <- test$p.value
      
      value <- aggregate(dat,by = list(group), mean)
      logFC <- value[1,2]-value[2,2]
      contrast <- paste0(value[1, 1], "-", value[2, 1])
      
      res <- data.frame(
        Symbol = symbol,
        logFC = logFC,
        contrast = contrast,
        pvalue = pvalue
      )
      return(res)
    })
    DEG <- bind_rows(DEG)
    rownames(DEG) <- NULL
    DEG <- na.omit(DEG)
    DEG$FDR <- p.adjust(DEG$pvalue, method = "BH")
    DEG$Sig <- fun_sig(DEG, logFC = 1, FDR = .05)
    return(DEG)
  })
  names(res_diff) <- map_chr(con, function(x)paste0(x, "-", "Other"))
  map(res_diff, ~table(.x$Sig))
  save(res_diff, file = "res_diff_proteomic_vsOther.Rdata")
}

load("res_diff_proteomic_vsOther.Rdata")

for(aim in names(res_diff)){
  
  data_vol <- res_diff[[aim]]
  
  color <- structure(names = c("Not", "Up", "Down"),
                     c("grey", "#ce181e", "#007cc0"))
  
  for_label <- data_vol%>%
    filter(Sig != "Not")%>%
    group_by(sign(logFC))%>%
    arrange(desc(abs(logFC)))%>%
    arrange(FDR)%>%
    slice(1:5)
  
  figure_title <- str_replace(aim, "-Other$", " vs Others")
  
  ggplot(data_vol) +
    aes(y = -log10(FDR), x = logFC) +
    geom_point(aes(color = Sig), shape = "circle", size = 3) +
    scale_color_manual(values = color)+
    geom_hline(yintercept = -log10(0.05), color = "black", alpha = 0.6, lty = 2, size = 1.2)+
    geom_vline(xintercept = c(-1, 1), color = "black", alpha = 0.6, lty = 2, size = 1.2)+
    ggrepel::geom_text_repel(
      data = for_label, 
      aes(label = Symbol),
      color = "black",
      min.segment.length = 0,
      box.padding = unit(.3, "lines"),
      segment.colour="black")+
    labs(x = "logFC",
         y = "-log10(FDR)",
         fill = "")+
    ggtitle(paste0("Proteomics ", figure_title))+
    theme_bw()+
    theme(plot.title = element_text(hjust = .5),
          legend.title = element_blank(),
          legend.position = "none",
          axis.text = element_text())
  ggsave(paste0("volcano_protein_", aim,".pdf"), width = 3.5, height = 3.5)
  
}
















