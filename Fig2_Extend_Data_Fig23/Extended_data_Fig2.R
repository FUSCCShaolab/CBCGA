rm(list = ls())

#library(dplyr)
#library(stringr)
#library(tidyr)
#library(data.table)
#library(survtype)
#library(ggplot2)
library(ggridges)
#library(escape)
library(maftools)
library(trackViewer)
library(export)
#library(BradleyTerry2)
library(tidyverse)

# Data preparation --------------------------------------------------------
load("/data/CBCGA.Extended_MergedData_V2.5_220722.Rdata")
CBCGAClin <- CBCGA_Cohort.Info
#rm(list = ls()[-c(20:21)])
rm(list = ls()[! ls() %in% c("CBCGAClin","CBCGA_WES_Somatic")])
load("/data/Fig2_S2_data.RData")

CBCGAClin <- janitor::clean_names(CBCGAClin)
colnames(CBCGAClin)[1] <- "Tumor_Sample_Barcode"
CBCGAClin <- subset(CBCGAClin, histological_type == "IDC")
CBCGAClin_WES <- subset(CBCGAClin, exome_sequencing_paired == "Yes")

TCGAClin$Tumor_Sample_Barcode <- str_sub(TCGAClin$sampleID, 1, 12)
TCGAClin_white <- subset(TCGAClin, Race.Category == "WHITE" & histological_type == "Infiltrating Ductal Carcinoma")
TCGAClin_WES <- subset(TCGAClin, Tumor_Sample_Barcode %in% TCGAMaf_df$Tumor_Sample_Barcode)
TCGAClin_WES_white <- subset(TCGAClin_white, Tumor_Sample_Barcode %in% TCGAMaf_df$Tumor_Sample_Barcode)

CBCGAClin_LumA <- subset(CBCGAClin_WES, intrinsic_subtype_pam50 == "LumA")
CBCGAMaf_df_LumA <- subset(CBCGA_WES_Somatic, Tumor_Sample_Barcode %in% CBCGAClin_LumA$Tumor_Sample_Barcode)
CBCGAMaf_LumA <- read.maf(CBCGAMaf_df_LumA)
TCGAClin_white_LumA <- subset(TCGAClin_WES_white, PAM50Call_RNAseq == "LumA")
TCGAMaf_df_LumA <- subset(TCGAMaf_df, Tumor_Sample_Barcode %in% TCGAClin_white_LumA$Tumor_Sample_Barcode)
TCGAMaf_LumA <- read.maf(TCGAMaf_df_LumA, vc_nonSyn = vc_nonSyn)

CBCGAvsTCGA_LumA <- coBarplot(TCGAMaf_LumA, CBCGAMaf_LumA,
                              genes = rev(genes),
                              m1Name = "TCGA Caucasian", m2Name = "CBCGA",
                              colors = mycol, yLims = c(50, 50), showPct = F
)
CBCGAvsTCGA_LumA_TCGA <- CBCGAvsTCGA_LumA$m1 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "TCGA")
CBCGAvsTCGA_LumA_CBCGA <- CBCGAvsTCGA_LumA$m2 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "CBCGA")
CBCGAvsTCGA_LumA_TCGA$Variant_Classification <- factor(CBCGAvsTCGA_LumA_TCGA$Variant_Classification,
                                                       levels = rev(c(
                                                         "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                         "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                       ))
)
CBCGAvsTCGA_LumA_CBCGA$Variant_Classification <- factor(CBCGAvsTCGA_LumA_CBCGA$Variant_Classification,
                                                        levels = rev(c(
                                                          "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                          "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                        ))
)
CBCGAvsTCGA_LumA_TCGA$name <- factor(CBCGAvsTCGA_LumA_TCGA$name, levels = rev(genes))
CBCGAvsTCGA_LumA_CBCGA$name <- factor(CBCGAvsTCGA_LumA_CBCGA$name, levels = rev(genes))
CBCGAvsTCGA_LumA_TCGA <- CBCGAvsTCGA_LumA_TCGA[order(CBCGAvsTCGA_LumA_TCGA$name, CBCGAvsTCGA_LumA_TCGA$Variant_Classification), ]
CBCGAvsTCGA_LumA_CBCGA <- CBCGAvsTCGA_LumA_CBCGA[order(CBCGAvsTCGA_LumA_CBCGA$name, CBCGAvsTCGA_LumA_CBCGA$Variant_Classification), ]
CBCGAvsTCGA_LumA_TCGA$value <- -CBCGAvsTCGA_LumA_TCGA$value


CBCGAClin_LumB <- subset(CBCGAClin_WES, intrinsic_subtype_pam50 == "LumB")
CBCGAMaf_df_LumB <- subset(CBCGA_WES_Somatic, Tumor_Sample_Barcode %in% CBCGAClin_LumB$Tumor_Sample_Barcode)
CBCGAMaf_LumB <- read.maf(CBCGAMaf_df_LumB)
TCGAClin_white_LumB <- subset(TCGAClin_WES_white, PAM50Call_RNAseq == "LumB")
TCGAMaf_df_LumB <- subset(TCGAMaf_df, Tumor_Sample_Barcode %in% TCGAClin_white_LumB$Tumor_Sample_Barcode)
TCGAMaf_LumB <- read.maf(TCGAMaf_df_LumB, vc_nonSyn = vc_nonSyn)

CBCGAvsTCGA_LumB <- coBarplot(TCGAMaf_LumB, CBCGAMaf_LumB,
                              genes = rev(genes),
                              m1Name = "TCGA Caucasian", m2Name = "CBCGA",
                              colors = mycol, yLims = c(50, 50), showPct = F
)
CBCGAvsTCGA_LumB_TCGA <- CBCGAvsTCGA_LumB$m1 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "TCGA")
CBCGAvsTCGA_LumB_CBCGA <- CBCGAvsTCGA_LumB$m2 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "CBCGA")
CBCGAvsTCGA_LumB_TCGA$Variant_Classification <- factor(CBCGAvsTCGA_LumB_TCGA$Variant_Classification,
                                                       levels = rev(c(
                                                         "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                         "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                       ))
)
CBCGAvsTCGA_LumB_CBCGA$Variant_Classification <- factor(CBCGAvsTCGA_LumB_CBCGA$Variant_Classification,
                                                        levels = rev(c(
                                                          "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                          "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                        ))
)
CBCGAvsTCGA_LumB_TCGA$name <- factor(CBCGAvsTCGA_LumB_TCGA$name, levels = rev(genes))
CBCGAvsTCGA_LumB_CBCGA$name <- factor(CBCGAvsTCGA_LumB_CBCGA$name, levels = rev(genes))
CBCGAvsTCGA_LumB_TCGA <- CBCGAvsTCGA_LumB_TCGA[order(CBCGAvsTCGA_LumB_TCGA$name, CBCGAvsTCGA_LumB_TCGA$Variant_Classification), ]
CBCGAvsTCGA_LumB_CBCGA <- CBCGAvsTCGA_LumB_CBCGA[order(CBCGAvsTCGA_LumB_CBCGA$name, CBCGAvsTCGA_LumB_CBCGA$Variant_Classification), ]
CBCGAvsTCGA_LumB_TCGA$value <- -CBCGAvsTCGA_LumB_TCGA$value


CBCGAClin_Her2 <- subset(CBCGAClin_WES, intrinsic_subtype_pam50 == "Her2")
CBCGAMaf_df_Her2 <- subset(CBCGA_WES_Somatic, Tumor_Sample_Barcode %in% CBCGAClin_Her2$Tumor_Sample_Barcode)
CBCGAMaf_Her2 <- read.maf(CBCGAMaf_df_Her2)
TCGAClin_white_Her2 <- subset(TCGAClin_WES_white, PAM50Call_RNAseq == "Her2")
TCGAMaf_df_Her2 <- subset(TCGAMaf_df, Tumor_Sample_Barcode %in% TCGAClin_white_Her2$Tumor_Sample_Barcode)
TCGAMaf_Her2 <- read.maf(TCGAMaf_df_Her2, vc_nonSyn = vc_nonSyn)

CBCGAvsTCGA_Her2 <- coBarplot(TCGAMaf_Her2, CBCGAMaf_Her2,
                              genes = rev(genes),
                              m1Name = "TCGA Caucasian", m2Name = "CBCGA",
                              colors = mycol, yLims = c(50, 50), showPct = F
)
CBCGAvsTCGA_Her2_TCGA <- CBCGAvsTCGA_Her2$m1 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "TCGA")
CBCGAvsTCGA_Her2_CBCGA <- CBCGAvsTCGA_Her2$m2 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "CBCGA")
CBCGAvsTCGA_Her2_TCGA$Variant_Classification <- factor(CBCGAvsTCGA_Her2_TCGA$Variant_Classification,
                                                       levels = rev(c(
                                                         "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                         "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                       ))
)
CBCGAvsTCGA_Her2_CBCGA$Variant_Classification <- factor(CBCGAvsTCGA_Her2_CBCGA$Variant_Classification,
                                                        levels = rev(c(
                                                          "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                          "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                        ))
)
CBCGAvsTCGA_Her2_TCGA$name <- factor(CBCGAvsTCGA_Her2_TCGA$name, levels = rev(genes))
CBCGAvsTCGA_Her2_CBCGA$name <- factor(CBCGAvsTCGA_Her2_CBCGA$name, levels = rev(genes))
CBCGAvsTCGA_Her2_TCGA <- CBCGAvsTCGA_Her2_TCGA[order(CBCGAvsTCGA_Her2_TCGA$name, CBCGAvsTCGA_Her2_TCGA$Variant_Classification), ]
CBCGAvsTCGA_Her2_CBCGA <- CBCGAvsTCGA_Her2_CBCGA[order(CBCGAvsTCGA_Her2_CBCGA$name, CBCGAvsTCGA_Her2_CBCGA$Variant_Classification), ]
CBCGAvsTCGA_Her2_TCGA$value <- -CBCGAvsTCGA_Her2_TCGA$value


CBCGAClin_Basal <- subset(CBCGAClin_WES, intrinsic_subtype_pam50 == "Basal")
CBCGAMaf_df_Basal <- subset(CBCGA_WES_Somatic, Tumor_Sample_Barcode %in% CBCGAClin_Basal$Tumor_Sample_Barcode)
CBCGAMaf_Basal <- read.maf(CBCGAMaf_df_Basal)
TCGAClin_white_Basal <- subset(TCGAClin_WES_white, PAM50Call_RNAseq == "Basal")
TCGAMaf_df_Basal <- subset(TCGAMaf_df, Tumor_Sample_Barcode %in% TCGAClin_white_Basal$Tumor_Sample_Barcode)
TCGAMaf_Basal <- read.maf(TCGAMaf_df_Basal, vc_nonSyn = vc_nonSyn)

CBCGAvsTCGA_Basal <- coBarplot(TCGAMaf_Basal, CBCGAMaf_Basal,
                               genes = rev(genes),
                               m1Name = "TCGA Caucasian", m2Name = "CBCGA",
                               colors = mycol, yLims = c(50, 50), showPct = F
)
CBCGAvsTCGA_Basal_TCGA <- CBCGAvsTCGA_Basal$m1 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "TCGA")
CBCGAvsTCGA_Basal_CBCGA <- CBCGAvsTCGA_Basal$m2 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "CBCGA")
CBCGAvsTCGA_Basal_TCGA$Variant_Classification <- factor(CBCGAvsTCGA_Basal_TCGA$Variant_Classification,
                                                        levels = rev(c(
                                                          "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                          "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                        ))
)
CBCGAvsTCGA_Basal_CBCGA$Variant_Classification <- factor(CBCGAvsTCGA_Basal_CBCGA$Variant_Classification,
                                                         levels = rev(c(
                                                           "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                           "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                         ))
)
CBCGAvsTCGA_Basal_TCGA$name <- factor(CBCGAvsTCGA_Basal_TCGA$name, levels = rev(genes))
CBCGAvsTCGA_Basal_CBCGA$name <- factor(CBCGAvsTCGA_Basal_CBCGA$name, levels = rev(genes))
CBCGAvsTCGA_Basal_TCGA <- CBCGAvsTCGA_Basal_TCGA[order(CBCGAvsTCGA_Basal_TCGA$name, CBCGAvsTCGA_Basal_TCGA$Variant_Classification), ]
CBCGAvsTCGA_Basal_CBCGA <- CBCGAvsTCGA_Basal_CBCGA[order(CBCGAvsTCGA_Basal_CBCGA$name, CBCGAvsTCGA_Basal_CBCGA$Variant_Classification), ]
CBCGAvsTCGA_Basal_TCGA$value <- -CBCGAvsTCGA_Basal_TCGA$value


CBCGAClin_Normal <- subset(CBCGAClin_WES, intrinsic_subtype_pam50 == "Normal")
CBCGAMaf_df_Normal <- subset(CBCGA_WES_Somatic, Tumor_Sample_Barcode %in% CBCGAClin_Normal$Tumor_Sample_Barcode)
CBCGAMaf_Normal <- read.maf(CBCGAMaf_df_Normal)
TCGAClin_white_Normal <- subset(TCGAClin_WES_white, PAM50Call_RNAseq == "Normal")
TCGAMaf_df_Normal <- subset(TCGAMaf_df, Tumor_Sample_Barcode %in% TCGAClin_white_Normal$Tumor_Sample_Barcode)
TCGAMaf_Normal <- read.maf(TCGAMaf_df_Normal, vc_nonSyn = vc_nonSyn)

CBCGAvsTCGA_Normal <- coBarplot(TCGAMaf_Normal, CBCGAMaf_Normal,
                                genes = rev(genes),
                                m1Name = "TCGA Caucasian", m2Name = "CBCGA",
                                colors = mycol, yLims = c(50, 50), showPct = F
)
CBCGAvsTCGA_Normal_TCGA <- CBCGAvsTCGA_Normal$m1 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "TCGA")
CBCGAvsTCGA_Normal_CBCGA <- CBCGAvsTCGA_Normal$m2 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "CBCGA")
CBCGAvsTCGA_Normal_TCGA$Variant_Classification <- factor(CBCGAvsTCGA_Normal_TCGA$Variant_Classification,
                                                         levels = rev(c(
                                                           "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                           "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                         ))
)
CBCGAvsTCGA_Normal_CBCGA$Variant_Classification <- factor(CBCGAvsTCGA_Normal_CBCGA$Variant_Classification,
                                                          levels = rev(c(
                                                            "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                            "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                          ))
)
CBCGAvsTCGA_Normal_TCGA$name <- factor(CBCGAvsTCGA_Normal_TCGA$name, levels = rev(genes))
CBCGAvsTCGA_Normal_CBCGA$name <- factor(CBCGAvsTCGA_Normal_CBCGA$name, levels = rev(genes))
CBCGAvsTCGA_Normal_TCGA <- CBCGAvsTCGA_Normal_TCGA[order(CBCGAvsTCGA_Normal_TCGA$name, CBCGAvsTCGA_Normal_TCGA$Variant_Classification), ]
CBCGAvsTCGA_Normal_CBCGA <- CBCGAvsTCGA_Normal_CBCGA[order(CBCGAvsTCGA_Normal_CBCGA$name, CBCGAvsTCGA_Normal_CBCGA$Variant_Classification), ]
CBCGAvsTCGA_Normal_TCGA$value <- -CBCGAvsTCGA_Normal_TCGA$value


CBCGAvsTCGA_AKT1_LumA <- rbind(CBCGAvsTCGA_LumA_CBCGA, CBCGAvsTCGA_LumA_TCGA) %>%
  filter(name == "AKT1") %>%
  mutate(Subtype = "LumA")
CBCGAvsTCGA_AKT1_LumB <- rbind(CBCGAvsTCGA_LumB_CBCGA, CBCGAvsTCGA_LumB_TCGA) %>%
  filter(name == "AKT1") %>%
  mutate(Subtype = "LumB")
CBCGAvsTCGA_AKT1_Her2 <- rbind(CBCGAvsTCGA_Her2_CBCGA, CBCGAvsTCGA_Her2_TCGA) %>%
  filter(name == "AKT1") %>%
  mutate(Subtype = "Her2")
CBCGAvsTCGA_AKT1_Basal <- rbind(CBCGAvsTCGA_Basal_CBCGA, CBCGAvsTCGA_Basal_TCGA) %>%
  filter(name == "AKT1") %>%
  mutate(Subtype = "Basal")
CBCGAvsTCGA_AKT1_Normal <- rbind(CBCGAvsTCGA_Normal_CBCGA, CBCGAvsTCGA_Normal_TCGA) %>%
  filter(name == "AKT1") %>%
  mutate(Subtype = "Normal")
CBCGAvsTCGA_AKT1 <- rbind(
  CBCGAvsTCGA_AKT1_LumA, CBCGAvsTCGA_AKT1_LumB, CBCGAvsTCGA_AKT1_Her2,
  CBCGAvsTCGA_AKT1_Basal, CBCGAvsTCGA_AKT1_Normal
)
CBCGAvsTCGA_AKT1$Subtype <- factor(CBCGAvsTCGA_AKT1$Subtype, levels = c("LumA", "LumB", "Her2", "Basal", "Normal"))


# Figure S2A --------------------------------------------------------------
pS2C <- ggplot() +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_LumA_CBCGA, position = "stack", stat = "identity"
  ) +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_LumA_TCGA, position = "stack", stat = "identity"
  ) +
  geom_hline(yintercept = 0) +
  guides(
    x = guide_axis(title = "Prevalence"),
    y = guide_axis(title = "Hugo symbol"),
    fill = guide_legend(title = "")
  ) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(labels = abs) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 10, colour = "black")
  )
graph2pdf(pS2C,file = '/results/FigS2A.pdf', width = 6, height = 6.9)

fisher_table <- data.frame()
for (i in genes) {
  tmp_CBCGA_mut <- genesToBarcodes(CBCGAMaf_LumA, genes = i, justNames = T)[[1]] %>% length()
  tmp_CBCGA_wt <- length(unique(CBCGAMaf_df_LumA$Tumor_Sample_Barcode)) - tmp_CBCGA_mut
  tmp_CBCGA_freq <- tmp_CBCGA_mut / (tmp_CBCGA_mut + tmp_CBCGA_wt)
  tmp_TCGA_mut <- genesToBarcodes(TCGAMaf_LumA, genes = i, justNames = T)[[1]] %>% length()
  tmp_TCGA_wt <- length(unique(TCGAMaf_df_LumA$Tumor_Sample_Barcode)) - tmp_TCGA_mut
  tmp_TCGA_freq <- tmp_TCGA_mut / (tmp_TCGA_mut + tmp_TCGA_wt)
  
  tmp_fisher <- fisher.test(matrix(c(tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_TCGA_mut, tmp_TCGA_wt), byrow = F, ncol = 2))
  tmp_pval <- tmp_fisher$p.value
  tmp_df <- c(i, tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_CBCGA_freq, tmp_TCGA_mut, tmp_TCGA_wt, tmp_TCGA_freq, tmp_pval)
  fisher_table <- rbind(fisher_table, tmp_df)
}
colnames(fisher_table) <- c("Symbol", "CBCGA_mut", "CBCGA_wt", "CBCGA_freq", "TCGA_mut", "TCGA_wt", "TCGA_freq", "pval")
fisher_table[, 2:8] <- fisher_table[, 2:8] %>% mutate_if(is.character, as.numeric)
fisher_table$FDR <- p.adjust(fisher_table$pval, method = "fdr")
# write.csv(fisher_table, 'Figure/Luminal A.csv', row.names = F)


# Figure S2B ---------------------------------------------------------------
pS2D <- ggplot() +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_LumB_CBCGA, position = "stack", stat = "identity"
  ) +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_LumB_TCGA, position = "stack", stat = "identity"
  ) +
  geom_hline(yintercept = 0) +
  guides(
    x = guide_axis(title = "Prevalence"),
    y = guide_axis(title = "Hugo symbol"),
    fill = guide_legend(title = "")
  ) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(labels = abs) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 10, colour = "black")
  )
graph2pdf(pS2D,file = '/results/FigS2B.pdf', width = 6, height = 6.9)


fisher_table <- data.frame()
for (i in genes) {
  tmp_CBCGA_mut <- genesToBarcodes(CBCGAMaf_LumB, genes = i, justNames = T)[[1]] %>% length()
  tmp_CBCGA_wt <- length(unique(CBCGAMaf_df_LumB$Tumor_Sample_Barcode)) - tmp_CBCGA_mut
  tmp_CBCGA_freq <- tmp_CBCGA_mut / (tmp_CBCGA_mut + tmp_CBCGA_wt)
  tmp_TCGA_mut <- genesToBarcodes(TCGAMaf_LumB, genes = i, justNames = T)[[1]] %>% length()
  tmp_TCGA_wt <- length(unique(TCGAMaf_df_LumB$Tumor_Sample_Barcode)) - tmp_TCGA_mut
  tmp_TCGA_freq <- tmp_TCGA_mut / (tmp_TCGA_mut + tmp_TCGA_wt)
  
  tmp_fisher <- fisher.test(matrix(c(tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_TCGA_mut, tmp_TCGA_wt), byrow = F, ncol = 2))
  tmp_pval <- tmp_fisher$p.value
  tmp_df <- c(i, tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_CBCGA_freq, tmp_TCGA_mut, tmp_TCGA_wt, tmp_TCGA_freq, tmp_pval)
  fisher_table <- rbind(fisher_table, tmp_df)
}
colnames(fisher_table) <- c("Symbol", "CBCGA_mut", "CBCGA_wt", "CBCGA_freq", "TCGA_mut", "TCGA_wt", "TCGA_freq", "pval")
fisher_table[, 2:8] <- fisher_table[, 2:8] %>% mutate_if(is.character, as.numeric)
fisher_table$FDR <- p.adjust(fisher_table$pval, method = "fdr")
# write.csv(fisher_table, 'Figure/Luminal B.csv', row.names = F)



# Figure S2C ---------------------------------------------------------------
pS2E <- ggplot() +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_Her2_CBCGA, position = "stack", stat = "identity"
  ) +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_Her2_TCGA, position = "stack", stat = "identity"
  ) +
  geom_hline(yintercept = 0) +
  guides(
    x = guide_axis(title = "Prevalence"),
    y = guide_axis(title = "Hugo symbol"),
    fill = guide_legend(title = "")
  ) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(limits = c(-80, 80), breaks = seq(-80, 80, 40), labels = abs) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 10, colour = "black")
  )
graph2pdf(pS2E,file = '/results/FigS2C.pdf', width = 6, height = 6.9)


fisher_table <- data.frame()
for (i in genes) {
  tmp_CBCGA_mut <- genesToBarcodes(CBCGAMaf_Her2, genes = i, justNames = T)[[1]] %>% length()
  tmp_CBCGA_wt <- length(unique(CBCGAMaf_df_Her2$Tumor_Sample_Barcode)) - tmp_CBCGA_mut
  tmp_CBCGA_freq <- tmp_CBCGA_mut / (tmp_CBCGA_mut + tmp_CBCGA_wt)
  tmp_TCGA_mut <- genesToBarcodes(TCGAMaf_Her2, genes = i, justNames = T)[[1]] %>% length()
  tmp_TCGA_wt <- length(unique(TCGAMaf_df_Her2$Tumor_Sample_Barcode)) - tmp_TCGA_mut
  tmp_TCGA_freq <- tmp_TCGA_mut / (tmp_TCGA_mut + tmp_TCGA_wt)
  
  tmp_fisher <- fisher.test(matrix(c(tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_TCGA_mut, tmp_TCGA_wt), byrow = F, ncol = 2))
  tmp_pval <- tmp_fisher$p.value
  tmp_df <- c(i, tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_CBCGA_freq, tmp_TCGA_mut, tmp_TCGA_wt, tmp_TCGA_freq, tmp_pval)
  fisher_table <- rbind(fisher_table, tmp_df)
}
colnames(fisher_table) <- c("Symbol", "CBCGA_mut", "CBCGA_wt", "CBCGA_freq", "TCGA_mut", "TCGA_wt", "TCGA_freq", "pval")
fisher_table[, 2:8] <- fisher_table[, 2:8] %>% mutate_if(is.character, as.numeric)
fisher_table$FDR <- p.adjust(fisher_table$pval, method = "fdr")
# write.csv(fisher_table, 'Figure/HER2 enriched.csv', row.names = F)





# Figure S2D ---------------------------------------------------------------
pS2F <- ggplot() +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_Basal_CBCGA, position = "stack", stat = "identity"
  ) +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_Basal_TCGA, position = "stack", stat = "identity"
  ) +
  geom_hline(yintercept = 0) +
  guides(
    x = guide_axis(title = "Prevalence"),
    y = guide_axis(title = "Hugo symbol"),
    fill = guide_legend(title = "")
  ) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(limits = c(-80, 80), breaks = seq(-80, 80, 40), labels = abs) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 10, colour = "black")
  )
graph2pdf(pS2F,file = '/results/FigS2D.pdf', width = 6, height = 6.9)


fisher_table <- data.frame()
for (i in genes) {
  tmp_CBCGA_mut <- genesToBarcodes(CBCGAMaf_Basal, genes = i, justNames = T)[[1]] %>% length()
  tmp_CBCGA_wt <- length(unique(CBCGAMaf_df_Basal$Tumor_Sample_Barcode)) - tmp_CBCGA_mut
  tmp_CBCGA_freq <- tmp_CBCGA_mut / (tmp_CBCGA_mut + tmp_CBCGA_wt)
  tmp_TCGA_mut <- genesToBarcodes(TCGAMaf_Basal, genes = i, justNames = T)[[1]] %>% length()
  tmp_TCGA_wt <- length(unique(TCGAMaf_df_Basal$Tumor_Sample_Barcode)) - tmp_TCGA_mut
  tmp_TCGA_freq <- tmp_TCGA_mut / (tmp_TCGA_mut + tmp_TCGA_wt)
  
  tmp_fisher <- fisher.test(matrix(c(tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_TCGA_mut, tmp_TCGA_wt), byrow = F, ncol = 2))
  tmp_pval <- tmp_fisher$p.value
  tmp_df <- c(i, tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_CBCGA_freq, tmp_TCGA_mut, tmp_TCGA_wt, tmp_TCGA_freq, tmp_pval)
  fisher_table <- rbind(fisher_table, tmp_df)
}
colnames(fisher_table) <- c("Symbol", "CBCGA_mut", "CBCGA_wt", "CBCGA_freq", "TCGA_mut", "TCGA_wt", "TCGA_freq", "pval")
fisher_table[, 2:8] <- fisher_table[, 2:8] %>% mutate_if(is.character, as.numeric)
fisher_table$FDR <- p.adjust(fisher_table$pval, method = "fdr")
# write.csv(fisher_table, 'Figure/Basal.csv', row.names = F)





# Figure S2E ---------------------------------------------------------------
pS2G <- ggplot() +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_Normal_CBCGA, position = "stack", stat = "identity"
  ) +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_Normal_TCGA, position = "stack", stat = "identity"
  ) +
  geom_hline(yintercept = 0) +
  guides(
    x = guide_axis(title = "Prevalence"),
    y = guide_axis(title = "Hugo symbol"),
    fill = guide_legend(title = "")
  ) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(limits = c(-60, 60), breaks = seq(-60, 60, 30), labels = abs) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 10, colour = "black")
  )
pS2G
# ggarrange(pS2C, pS2D, pS2E, pS2F, pS2G,
#           ncol = 2, nrow = 3, labels = c("C", "D", "E", "F", "G"),
#           legend = "none"
# )
graph2pdf(pS2G,file = '/results/FigS2E.pdf', width = 6, height = 6.9)


fisher_table <- data.frame()
for (i in genes) {
  tmp_CBCGA_mut <- genesToBarcodes(CBCGAMaf_Normal, genes = i, justNames = T)[[1]] %>% length()
  tmp_CBCGA_wt <- length(unique(CBCGAMaf_df_Normal$Tumor_Sample_Barcode)) - tmp_CBCGA_mut
  tmp_CBCGA_freq <- tmp_CBCGA_mut / (tmp_CBCGA_mut + tmp_CBCGA_wt)
  tmp_TCGA_mut <- genesToBarcodes(TCGAMaf_Normal, genes = i, justNames = T)[[1]] %>% length()
  tmp_TCGA_wt <- length(unique(TCGAMaf_df_Normal$Tumor_Sample_Barcode)) - tmp_TCGA_mut
  tmp_TCGA_freq <- tmp_TCGA_mut / (tmp_TCGA_mut + tmp_TCGA_wt)
  
  tmp_fisher <- fisher.test(matrix(c(tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_TCGA_mut, tmp_TCGA_wt), byrow = F, ncol = 2))
  tmp_pval <- tmp_fisher$p.value
  tmp_df <- c(i, tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_CBCGA_freq, tmp_TCGA_mut, tmp_TCGA_wt, tmp_TCGA_freq, tmp_pval)
  fisher_table <- rbind(fisher_table, tmp_df)
}
colnames(fisher_table) <- c("Symbol", "CBCGA_mut", "CBCGA_wt", "CBCGA_freq", "TCGA_mut", "TCGA_wt", "TCGA_freq", "pval")
fisher_table[, 2:8] <- fisher_table[, 2:8] %>% mutate_if(is.character, as.numeric)
fisher_table$FDR <- p.adjust(fisher_table$pval, method = "fdr")
# write.csv(fisher_table, 'Figure/Normal.csv', row.names = F)
