rm(list = ls())

library(tidyverse)
library(clusterProfiler)
library(DescTools)
library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)
library(xlsx)
library(ggsci)

# 
load("CBCGA.Extended_MergedData_V2.6_220802.Rdata")
load("Data/meta_pathway.Rdata")

# 
info_subtype <- CBCGA_Cohort.Info[, c(1, 10:11, 44)]
info_subtype$`Intrinsic subtype (AIMS)`[info_subtype$`Intrinsic subtype (AIMS)` == "NA"] <- NA
# info_subtype <- info_subtype[!is.na(info_subtype$`Intrinsic subtype (PAM50)`), ]
names(info_subtype) <- c("PatientCode", "PAM50 subtypes", "AIMS subtypes", "IHC subtypes")

# 
if(TRUE){
  data <- CBCGA.Extended_lip
  info <- CBCGA_lip_anno
  info <- rownames_to_column(info, var = "Symbol")
  
  id_dat <- names(data)
  id_dat <- id_dat[str_detect(id_dat, "_T$")]
  table(duplicated(str_extract(id_dat, "[A-Z]*_")))
  
  info_subtype$datID <- id_dat[match(info_subtype$PatientCode, str_extract(id_dat, "[A-Z]{4,4}"))]
  info_subtype <- info_subtype[!is.na(info_subtype$datID), ]
  info_subtype <- info_subtype[!is.na(info_subtype$`PAM50 subtypes`), ]
  
  data <- data[info_subtype$datID]
  data <- data[match(info_subtype$datID, names(data))]
  identical(info_subtype$datID, names(data))
  
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
    con_nm <- paste(gg, collapse = "-")
    con_info <- info_subtype
    con_info$PAM50 <- ifelse(con_info$PAM50 == gg, gg, "Other")
    con_info$PAM50 <- factor(con_info$PAM50, levels = c(gg, "Other"))
    dat <- data[names(data)%in%con_info$datID]
    
    DEG <- map(seq_along(rownames(dat)), function(x){
      symbol <- rownames(dat)[x]
      dat <- unlist(dat[x, ])
      group <- con_info$PAM50
      test <- wilcox.test(dat~group)
      pvalue <- test$p.value
      
      value <- aggregate(dat, by = list(group), mean)
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
    DEG$Sig <- fun_sig(DEG, logFC = 1, FDR = 0.05)
    return(DEG)
    
  })
  names(res_diff) <- map_chr(con, function(x)paste0(x, "-", "Other"))
  map(res_diff, ~table(.x$Sig))
  save(res_diff, file = "res_diff_metaLipid_vsOther_PAM50.Rdata")
}

load("res_diff_metaLipid_vsOther_PAM50.Rdata")
for(aim in names(res_diff)){
  
  data_vol <- res_diff[[aim]]
  data_vol$Category <- info$Lipid.super.class[match(data_vol$Symbol, info$Symbol)]
  data_vol$group2 <- ifelse(data_vol$Sig == "Not", data_vol$Sig, data_vol$Category)
  data_vol$group3 <- ifelse(data_vol$Sig == "Not", data_vol$Sig, info$Subclass[match(data_vol$Symbol, info$Symbol)])
  

  color <- structure(c("#E82029", "#1F78B6", "#36A754", "#593D8E", "#E8E8E8"), 
                     names = c("Fatty acyls [FA]",  "Glycerolipids [GL]", "Glycerophospholipids [GP]", "Sphingolipids [SP]", "Not"))
  data_vol$color <- color[match(data_vol$group2, names(color))]
  
  figure_title <- str_replace(aim, "-Other$", " vs Others")
  
  ggplot(data_vol) +
    aes(y = -log10(FDR), x = logFC) +
    geom_point(aes(color = group2), shape = "circle", size = 3) +
    scale_color_manual(values = color)+
    geom_hline(yintercept = -log10(0.05), color = "black", alpha = 0.6, lty = 2, size = 1.2)+
    geom_vline(xintercept = c(-1, 1), color = "black", alpha = 0.6, lty = 2, size = 1.2)+
    labs(x = "logFC",
         y = "-log10(FDR)",
         fill = "")+
    ggtitle(paste0("Metabolic Lipids ", figure_title))+
    theme_bw()+
    theme(plot.title = element_text(hjust = .5),
          legend.title = element_blank(),
          legend.position = "none",
          axis.text = element_text())
  ggsave(paste0("volcano_meta_lipid", aim,".pdf"), 
         width = 3.5, height = 3.5)
  
}

