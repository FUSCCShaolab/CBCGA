library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(data.table)
library(survtype)
library(ggplot2)
library(ggridges)
library(escape)
library(maftools)
library(trackViewer)
library(export)
library(BradleyTerry2)



# Data preparation --------------------------------------------------------
load("E:/FDU_SCI/CBCGA/WES/R/Figure2/Data/CBCGA.Extended_MergedData_V2.5_220722.Rdata")
CBCGAClin <- CBCGA_Cohort.Info
rm(list = ls()[-c(11, 20:21)])
load("E:/FDU_SCI/CBCGA/WES/R/Figure2/Data/Other data for Figure 2.RData")

CBCGAClin <- janitor::clean_names(CBCGAClin)
colnames(CBCGAClin)[1] <- "Tumor_Sample_Barcode"
CBCGAClin <- subset(CBCGAClin, histological_type == "IDC")
CBCGAClin_WES <- subset(CBCGAClin, exome_sequencing_paired == "Yes")

TCGAClin$Tumor_Sample_Barcode <- str_sub(TCGAClin$sampleID, 1, 12)
TCGAClin_white <- subset(TCGAClin, Race.Category == "WHITE" & histological_type == "Infiltrating Ductal Carcinoma")
TCGAClin_WES <- subset(TCGAClin, Tumor_Sample_Barcode %in% TCGAMaf_df$Tumor_Sample_Barcode)
TCGAClin_WES_white <- subset(TCGAClin_white, Tumor_Sample_Barcode %in% TCGAMaf_df$Tumor_Sample_Barcode)





# Figure 2A ---------------------------------------------------------------
genes <- c(
  "TP53", "PIK3CA", "GATA3", "MAP3K1", "KMT2C", "AKT1", "PTEN", "NF1",
  "ARID1A", "SF3B1", "FAT4", "LRP1B", "MAP2K4"
)
mycol <- c(
  "Splice_Site" = "#6b4694", "Missense_Mutation" = "#4d93cd", "Nonsense_Mutation" = "#e63a3a",
  "Frame_Shift_Ins" = "#fab71b", "Frame_Shift_Del" = "#fab71b", "In_Frame_Ins" = "#ecd71e",
  "In_Frame_Del" = "#ecd71e", "Nonstop_Mutation" = "#ab5b9e", "Translation_Start_Site" = "#018b38",
  "Multi_Hit" = "#cecece"
)
vc_nonSyn <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del",
  "Splice_Site", "Splice_Region", "Nonstop_Mutation", "Translation_Start_Site"
)


CBCGAMaf_df <- subset(CBCGA_WES_Somatic, Tumor_Sample_Barcode %in% CBCGAClin$Tumor_Sample_Barcode)
CBCGAMaf <- read.maf(CBCGAMaf_df)
TCGAMaf_df_white <- subset(TCGAMaf_df, Tumor_Sample_Barcode %in% TCGAClin_white$Tumor_Sample_Barcode)
TCGAMaf_white <- read.maf(TCGAMaf_df_white, vc_nonSyn = vc_nonSyn)

CBCGAvsTCGA <- coBarplot(TCGAMaf_white, CBCGAMaf,
  genes = rev(genes),
  m1Name = "TCGA Caucasian", m2Name = "CBCGA",
  colors = mycol, yLims = c(50, 50), showPct = F
)
CBCGAvsTCGA_TCGA <- CBCGAvsTCGA$m1 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "TCGA")
CBCGAvsTCGA_CBCGA <- CBCGAvsTCGA$m2 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "CBCGA")
CBCGAvsTCGA_TCGA$Variant_Classification <- factor(CBCGAvsTCGA_TCGA$Variant_Classification,
  levels = rev(c(
    "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
    "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
  ))
)
CBCGAvsTCGA_CBCGA$Variant_Classification <- factor(CBCGAvsTCGA_CBCGA$Variant_Classification,
  levels = rev(c(
    "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
    "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
  ))
)
CBCGAvsTCGA_TCGA$name <- factor(CBCGAvsTCGA_TCGA$name, levels = rev(genes))
CBCGAvsTCGA_CBCGA$name <- factor(CBCGAvsTCGA_CBCGA$name, levels = rev(genes))
CBCGAvsTCGA_TCGA <- CBCGAvsTCGA_TCGA[order(CBCGAvsTCGA_TCGA$name, CBCGAvsTCGA_TCGA$Variant_Classification), ]
CBCGAvsTCGA_CBCGA <- CBCGAvsTCGA_CBCGA[order(CBCGAvsTCGA_CBCGA$name, CBCGAvsTCGA_CBCGA$Variant_Classification), ]
CBCGAvsTCGA_TCGA$value <- -CBCGAvsTCGA_TCGA$value

p2a <- ggplot() +
  geom_bar(aes(name, value, fill = Variant_Classification),
    data = CBCGAvsTCGA_CBCGA, position = "stack", stat = "identity"
  ) +
  geom_bar(aes(name, value, fill = Variant_Classification),
    data = CBCGAvsTCGA_TCGA, position = "stack", stat = "identity"
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
p2a





# Figure 2B ---------------------------------------------------------------
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


p2b <- ggplot() +
  geom_bar(aes(Subtype, value, fill = Variant_Classification),
    data = subset(CBCGAvsTCGA_AKT1, Cohort == "CBCGA"), position = "stack", stat = "identity"
  ) +
  geom_bar(aes(Subtype, value, fill = Variant_Classification),
    data = subset(CBCGAvsTCGA_AKT1, Cohort == "TCGA"), position = "stack", stat = "identity"
  ) +
  geom_hline(yintercept = 0) +
  guides(
    x = guide_axis(title = "", angle = 45),
    y = guide_axis(title = "Prevalence"),
    fill = guide_legend(title = "")
  ) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(labels = abs) +
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 10, colour = "black")
  )
p2b





# Figure 2c ---------------------------------------------------------------
CBCGAMaf_df_LumA3 <- CBCGAMaf_LumA@data
CBCGAMaf_df_LumA3$Chromosome <- gsub("chr", "", CBCGAMaf_df_LumA3$Chromosome)
CBCGAMaf_df_LumA3$mutation_id <- paste(CBCGAMaf_df_LumA3$Tumor_Sample_Barcode, CBCGAMaf_df_LumA3$Chromosome,
  CBCGAMaf_df_LumA3$Start_Position, CBCGAMaf_df_LumA3$Reference_Allele,
  sep = ":"
)


DriverGene <- c(
  "TP53", "AKT1", "PIK3CA", "MAP3K1", "PTEN", "MAP2K4", "CBFB", "GATA3", "KMT2C", "CDH1", "RUNX1", "FOXA1", "CDKN1B", "GPS2", "ARID1A",
  "ZFP36L1", "NF1", "PIK3R1", "SF3B1", "CTCF", "TBX3", "SMAD4"
)
CBCGA_pyclone_mutation_LumA <- subset(CBCGA_pyclone_mutation, Sample %in% CBCGAClin_LumA$Tumor_Sample_Barcode)
CBCGA_pyclone_mutation_LumA <- subset(CBCGA_pyclone_mutation_LumA, Gene.refGene %in% DriverGene)
CBCGA_pyclone_mutation_LumA <- CBCGA_pyclone_mutation_LumA[intersect(rownames(CBCGA_pyclone_mutation_LumA), unique(CBCGAMaf_df_LumA3$mutation_id)), ]
TCGA_pyclone <- fread("E:/FDU_SCI/CBCGA/Race/CBCGARaceR/MutationalTiming/Pyclone/ITH_mutations_final.csv")
TCGA_pyclone_mutation_LumA <- subset(TCGA_pyclone, Sample %in% TCGAClin_white_LumA$Tumor_Sample_Barcode & Gene.refGene %in% DriverGene)


genePair <- combn(DriverGene, 2) %>%
  t() %>%
  data.frame()
genePair$CBCGA <- NA
genePair$TCGA <- NA
if (T) {
  for (i in 1:nrow(genePair)) {
    tmp_ID1_CBCGA <- subset(CBCGA_pyclone_mutation_LumA, Gene.refGene == genePair[i, 1])$Sample %>% unique()
    tmp_ID2_CBCGA <- subset(CBCGA_pyclone_mutation_LumA, Gene.refGene == genePair[i, 2])$Sample %>% unique()
    tmp_ID_CBCGA <- intersect(tmp_ID1_CBCGA, tmp_ID2_CBCGA) %>% length()
    genePair[i, "CBCGA"] <- tmp_ID_CBCGA

    tmp_ID1_TCGA <- subset(TCGA_pyclone_mutation_LumA, Gene.refGene == genePair[i, 1])$Sample %>% unique()
    tmp_ID2_TCGA <- subset(TCGA_pyclone_mutation_LumA, Gene.refGene == genePair[i, 2])$Sample %>% unique()
    tmp_ID_TCGA <- intersect(tmp_ID1_TCGA, tmp_ID2_TCGA) %>% length()
    genePair[i, "TCGA"] <- tmp_ID_TCGA
  }


  genePair_CBCGA <- subset(genePair, CBCGA >= 3)
  for (i in 1:nrow(genePair_CBCGA)) {
    tmp_ID1_CBCGA <- subset(CBCGA_pyclone_mutation_LumA, Gene.refGene == genePair_CBCGA[i, 1])$Sample %>% unique()
    tmp_ID2_CBCGA <- subset(CBCGA_pyclone_mutation_LumA, Gene.refGene == genePair_CBCGA[i, 2])$Sample %>% unique()
    tmp_ID_CBCGA <- intersect(tmp_ID1_CBCGA, tmp_ID2_CBCGA)

    tmp_win <- NULL
    for (j in tmp_ID_CBCGA) {
      tmp_VAF1_CBCGA <- subset(CBCGA_pyclone_mutation_LumA, Gene.refGene == genePair_CBCGA[i, 1] & Sample == j)$PreClusterCCF %>% mean()
      tmp_VAF2_CBCGA <- subset(CBCGA_pyclone_mutation_LumA, Gene.refGene == genePair_CBCGA[i, 2] & Sample == j)$PreClusterCCF %>% mean()
      tmp_win2 <- ifelse(tmp_VAF1_CBCGA - tmp_VAF2_CBCGA > 0.05, "CBCGA_win1",
        ifelse(tmp_VAF1_CBCGA - tmp_VAF2_CBCGA < -0.05, "CBCGA_win2", "CBCGA_amb")
      )
      tmp_win <- c(tmp_win, tmp_win2)
    }
    tmp_df_CBCGA <- table(tmp_win) %>% data.frame()
    if ("CBCGA_win1" %in% tmp_df_CBCGA$tmp_win) {
      genePair_CBCGA[i, "CBCGA_win1"] <- tmp_df_CBCGA[tmp_df_CBCGA$tmp_win == "CBCGA_win1", "Freq"]
    } else {
      genePair_CBCGA[i, "CBCGA_win1"] <- 0
    }
    if ("CBCGA_win2" %in% tmp_df_CBCGA$tmp_win) {
      genePair_CBCGA[i, "CBCGA_win2"] <- tmp_df_CBCGA[tmp_df_CBCGA$tmp_win == "CBCGA_win2", "Freq"]
    } else {
      genePair_CBCGA[i, "CBCGA_win2"] <- 0
    }
    if ("CBCGA_amb" %in% tmp_df_CBCGA$tmp_win) {
      genePair_CBCGA[i, "CBCGA_amb"] <- tmp_df_CBCGA[tmp_df_CBCGA$tmp_win == "CBCGA_amb", "Freq"]
    } else {
      genePair_CBCGA[i, "CBCGA_amb"] <- 0
    }
  }



  genePair_TCGA <- subset(genePair, TCGA >= 3)
  for (i in 1:nrow(genePair_TCGA)) {
    tmp_ID1_TCGA <- subset(TCGA_pyclone_mutation_LumA, Gene.refGene == genePair_TCGA[i, 1])$Sample %>% unique()
    tmp_ID2_TCGA <- subset(TCGA_pyclone_mutation_LumA, Gene.refGene == genePair_TCGA[i, 2])$Sample %>% unique()
    tmp_ID_TCGA <- intersect(tmp_ID1_TCGA, tmp_ID2_TCGA)

    tmp_win <- NULL
    for (j in tmp_ID_TCGA) {
      tmp_VAF1_TCGA <- subset(TCGA_pyclone_mutation_LumA, Gene.refGene == genePair_TCGA[i, 1] & Sample == j)$PreClusterCCF %>% mean()
      tmp_VAF2_TCGA <- subset(TCGA_pyclone_mutation_LumA, Gene.refGene == genePair_TCGA[i, 2] & Sample == j)$PreClusterCCF %>% mean()
      tmp_win2 <- ifelse(tmp_VAF1_TCGA - tmp_VAF2_TCGA > 0.05, "TCGA_win1",
        ifelse(tmp_VAF1_TCGA - tmp_VAF2_TCGA < -0.05, "TCGA_win2", "TCGA_amb")
      )
      tmp_win <- c(tmp_win, tmp_win2)
    }
    tmp_df_TCGA <- table(tmp_win) %>% data.frame()
    if ("TCGA_win1" %in% tmp_df_TCGA$tmp_win) {
      genePair_TCGA[i, "TCGA_win1"] <- tmp_df_TCGA[tmp_df_TCGA$tmp_win == "TCGA_win1", "Freq"]
    } else {
      genePair_TCGA[i, "TCGA_win1"] <- 0
    }
    if ("TCGA_win2" %in% tmp_df_TCGA$tmp_win) {
      genePair_TCGA[i, "TCGA_win2"] <- tmp_df_TCGA[tmp_df_TCGA$tmp_win == "TCGA_win2", "Freq"]
    } else {
      genePair_TCGA[i, "TCGA_win2"] <- 0
    }
    if ("TCGA_amb" %in% tmp_df_TCGA$tmp_win) {
      genePair_TCGA[i, "TCGA_amb"] <- tmp_df_TCGA[tmp_df_TCGA$tmp_win == "TCGA_amb", "Freq"]
    } else {
      genePair_TCGA[i, "TCGA_amb"] <- 0
    }
  }



  # Bradley-Terry2
  colnames(genePair_CBCGA)[c(1, 2, 5, 6)] <- c("Gene1", "Gene2", "Gene1_wins", "Gene2_wins")
  genePair_CBCGA <- genePair_CBCGA[, c(1, 2, 5, 6)]
  genePair_CBCGA$Gene1 <- factor(genePair_CBCGA$Gene1, levels = unique(c(genePair_CBCGA$Gene1, genePair_CBCGA$Gene2)))
  genePair_CBCGA$Gene2 <- factor(genePair_CBCGA$Gene2, levels = unique(c(as.character(genePair_CBCGA$Gene1), genePair_CBCGA$Gene2)))
  CBCGA_btData <- BTm(cbind(Gene1_wins, Gene2_wins), Gene1, Gene2, data = genePair_CBCGA)
  CBCGA_btData_ordering <- BTabilities(CBCGA_btData) %>%
    data.frame() %>%
    rownames_to_column("Symbol")
  CBCGA_btData_ordering <- subset(CBCGA_btData_ordering, s.e. <= 2)

  colnames(genePair_TCGA)[c(1, 2, 5, 6)] <- c("Gene1", "Gene2", "Gene1_wins", "Gene2_wins")
  genePair_TCGA <- genePair_TCGA[, c(1, 2, 5, 6)]
  genePair_TCGA$Gene1 <- factor(genePair_TCGA$Gene1, levels = unique(c(genePair_TCGA$Gene1, genePair_TCGA$Gene2)))
  genePair_TCGA$Gene2 <- factor(genePair_TCGA$Gene2, levels = unique(c(as.character(genePair_TCGA$Gene1), genePair_TCGA$Gene2)))
  TCGA_btData <- BTm(cbind(Gene1_wins, Gene2_wins), Gene1, Gene2, data = genePair_TCGA)
  TCGA_btData_ordering <- BTabilities(TCGA_btData) %>%
    data.frame() %>%
    rownames_to_column("Symbol")
  TCGA_btData_ordering <- subset(TCGA_btData_ordering, s.e. <= 2)

  CBCGA_btData_ordering$Pathway <- c(
    "TP53", "PI3K", "PI3K", "Other", "PI3K", "Transcription factor", "Transcription factor",
    "Epigenetic modifiers", "Transcription factor", "Cell cycle", "Splicing"
  )
  TCGA_btData_ordering$Pathway <- c(
    "TP53", "PI3K", "Transcription factor", "Transcription factor", "Other", "PI3K", "Other",
    "Transcription factor", "RTK", "Epigenetic modifiers"
  )

  mycol2 <- c("#A94543", "#1F97C6", "#75A992", "#D68A46", "#334D54", "#666393", "#7F796B", "#F9BD1D")
  names(mycol2) <- c(
    "PI3K", "Cell cycle", "Transcription factor", "Epigenetic modifiers",
    "Splicing", "Other", "RTK", "TP53"
  )
  p2g1 <- ggplot(
    CBCGA_btData_ordering,
    aes(reorder(factor(Symbol), ability),
      y = ability,
      ymin = (ability - s.e.), ymax = (ability + s.e.)
    )
  ) +
    geom_pointrange(
      size = 0.75,
      aes(color = Pathway)
    ) +
    scale_color_manual(values = mycol2) +
    scale_y_reverse() +
    guides(
      x = guide_axis(title = "Point Estimate + 95% CI"),
      y = guide_axis(title = ""),
      color = guide_legend(title = "")
    ) +
    coord_flip() +
    theme_classic() +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text = element_text(colour = "black", size = 10),
      axis.title = element_text(size = 12)
    )
}
p2g1


CBCGA_btData_ordering <- CBCGA_btData_ordering[order(CBCGA_btData_ordering$ability, decreasing = T), ]
CBCGA_pyclone_mutation_LumA2 <- subset(CBCGA_pyclone_mutation_LumA, Gene.refGene %in% CBCGA_btData_ordering$Symbol)
CBCGA_pyclone_mutation_LumA2$Gene.refGene <- factor(CBCGA_pyclone_mutation_LumA2$Gene.refGene, levels = rev(CBCGA_btData_ordering$Symbol))
# CBCGA_pyclone_mutation_LumA2 <- merge(CBCGA_pyclone_mutation_LumA2, CBCGA_btData_ordering[, 5, drop = F], by.x = "Gene.refGene", by.y = "Symbol", all.x = T)
p2g2 <- ggplot(CBCGA_pyclone_mutation_LumA2, aes(PreClusterCCF, Gene.refGene)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01) +
  scale_fill_manual(values = mycol2) +
  guides(
    y = guide_axis(title = ""),
    fill = guide_legend(title = "")
  ) +
  theme_classic()

ggarrange(p2g2, p2g1, ncol = 2, nrow = 1, common.legend = T, legend = "right")




# Figure 2D ---------------------------------------------------------------

TCGAClin_IDC <- TCGAClin %>%
  filter(
    Race.Category == "WHITE",
    histological_type == "Infiltrating Ductal Carcinoma",
    PAM50Call_RNAseq != ""
  )
CBCGAClin_IDC <- CBCGA_Cohort.Info %>%
  clean_names() %>%
  filter(
    histological_type == "IDC",
    is.na(intrinsic_subtype_pam50) == F
  )
tmp_df <- rbind(
  data.frame(clinical_subtype = CBCGAClin_IDC$clinical_subtype, pam50 = CBCGAClin_IDC$intrinsic_subtype_pam50, cohort = "cbcga"),
  data.frame(clinical_subtype = TCGAClin_IDC$Clinical_Subtype, pam50 = TCGAClin_IDC$PAM50Call_RNAseq, cohort = "tcga")
)

tmp_df2 <- rbind(tmp_df %>% mutate(clinical_subtype = "All"), tmp_df)
table(tmp_df2$pam50, tmp_df2$clinical_subtype, tmp_df2$cohort) %>%
  data.frame() %>%
  filter(Var2 != "Unknown") %>%
  mutate(
    Var1 = factor(Var1, levels = c("LumA", "LumB", "Her2", "Basal", "Normal")),
    Var2 = factor(Var2, levels = c("All", "HR+HER2-", "HR+HER2+", "HR-HER2+", "TNBC"))
  ) %>%
  ggplot(aes(Var3, Freq, fill = Var1)) +
  geom_bar(position = "fill", stat = "identity", color = "black") +
  scale_y_continuous("Proportion", expand = c(0, 0)) +
  scale_fill_manual(values = c("#1D76BB", "#76CFE6", "#6D58A6", "#E01D2D", "#CDCFD0")) +
  facet_wrap(~Var2, ncol = 5) +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 0.2358),
    axis.line.x = element_blank(),
    axis.ticks = element_line(linewidth = 0.2358, colour = "black"),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10, colour = "black"),
    axis.text.x = element_blank()
  )

# Figure 2E-F-----------------------------------------------------------------
for ( k in c("HRposHER2neg","HRposHER2pos","HRnegHER2pos","HRnegHER2neg")){
  
  CBCGA = subtype_CBCGA[[k]]
  TCGA = subtype_TCGA[[k]]
  
  tmp = intersect(cnvlist,rownames(TCGA_cnv))
  TCGA_cnv_amp = TCGA_cnv_loss = TCGA_cnv[tmp,TCGA]
  thre_amp = thre_loss = thre[tmp,CBCGA]
  
  for (i in rownames(TCGA_cnv_amp)){
    mat1_loss = mat1_amp = TCGA_cnv_amp[i,]
    mat1_amp[mat1_amp != 2 ] = 0 ; mat1_amp[mat1_amp == 2 ] = 1
    mat1_loss[mat1_loss != -2 ] = 0 ; mat1_loss[mat1_loss == -2 ] = 1
    
    TCGA_cnv_amp[i,] = mat1_amp
    TCGA_cnv_loss[i,] = mat1_loss
    
    mat2_loss = mat2_amp = thre_amp[i,]
    mat2_amp[mat2_amp != 2 ] = 0 ; mat2_amp[mat2_amp == 2 ] = 1
    mat2_loss[mat2_loss != -2 ] = 0 ; mat2_loss[mat2_loss == -2 ] = 1
    
    thre_amp[i,] = mat2_amp
    thre_loss[i,] = mat2_loss
    
    print(i)
  }
  
  #####compare
  CBCGAVStcga.loss.res = CBCGAVStcga.gain.res = as.data.frame( matrix(nrow = nrow(TCGA_cnv_amp), ncol = 6) )
  rownames(CBCGAVStcga.loss.res) = rownames(CBCGAVStcga.gain.res)  = row.names(TCGA_cnv_amp)
  colnames(CBCGAVStcga.loss.res) = colnames(CBCGAVStcga.gain.res) = c("CBCGA","TCGA","p.val","adj.p","CBCGA_case","TCGA_case")
  
  CBCGAVStcga.gain.res$CBCGA = apply(thre_amp,1,function(x){sum(x)/(length(x))}) 
  CBCGAVStcga.gain.res$CBCGA_case = apply(thre_amp,1,function(x){sum(x)}) 
  CBCGAVStcga.gain.res$TCGA = apply(TCGA_cnv_amp,1,function(x){sum(x)/(length(x))})
  CBCGAVStcga.gain.res$TCGA_case = apply(TCGA_cnv_amp,1,function(x){sum(x)})
  
  CBCGAVStcga.loss.res$CBCGA = apply(thre_loss,1,function(x){sum(x)/(length(x))})
  CBCGAVStcga.loss.res$CBCGA_case = apply(thre_loss,1,function(x){sum(x)})
  CBCGAVStcga.loss.res$TCGA = apply(TCGA_cnv_loss,1,function(x){sum(x)/(length(x))})
  CBCGAVStcga.loss.res$TCGA_case = apply(TCGA_cnv_loss,1,function(x){sum(x)})
  
  ## gain
  for ( i in rownames(CBCGAVStcga.gain.res)){
    tmp1 = t(rbind(rep("CBCGA",ncol(thre_amp)), thre_amp[i,]))
    tmp2= t(rbind(rep("TCGA",ncol(TCGA_cnv_amp)), TCGA_cnv_amp[i,]))
    tmp = as.data.frame(rbind(tmp1,tmp2))
    colnames(tmp) = c("cohort","gene")
    mytab = table(tmp$cohort,tmp$gene)
    if (dim(mytab)[2] != 1 ){
      #if( all(as.data.frame(mytab)[,3] > 5) ){p = chisq.test(mytab)} else {p = fisher.test(mytab)}
      p = fisher.test(mytab)
      CBCGAVStcga.gain.res[i,3] <- p$p.value
    }
  }
  CBCGAVStcga.gain.res[,4] = p.adjust(CBCGAVStcga.gain.res[,3],method = "fdr")
  CBCGAVStcga.gain.res = CBCGAVStcga.gain.res[!is.na(CBCGAVStcga.gain.res$p.val),]
  #write.csv(CBCGAVStcga.gain.res,file = paste0("/results/CBCGAVStcga_WHITE_gain_res_",k,".csv"))
  
  ## loss
  for ( i in rownames(CBCGAVStcga.loss.res)){
    tmp1 = t(rbind(rep("CBCGA",ncol(thre_loss)), thre_loss[i,]))
    tmp2= t(rbind(rep("TCGA",ncol(TCGA_cnv_loss)), TCGA_cnv_loss[i,]))
    tmp = as.data.frame(rbind(tmp1,tmp2))
    colnames(tmp) = c("cohort","gene")
    mytab = table(tmp$cohort,tmp$gene)
    if (dim(mytab)[2] != 1 ){
      #if( all(as.data.frame(mytab)[,3] > 5) ){p = chisq.test(mytab)} else {p = fisher.test(mytab)}
      p = fisher.test(mytab)
      CBCGAVStcga.loss.res[i,3] <- p$p.value
    }
  }
  CBCGAVStcga.loss.res[,4] = p.adjust(CBCGAVStcga.loss.res[,3],method = "fdr")
  CBCGAVStcga.loss.res = CBCGAVStcga.loss.res[!is.na(CBCGAVStcga.loss.res$p.val),]
  #write.csv(CBCGAVStcga.loss.res,file = paste0("/results/CBCGAVStcga_WHITE_loss_res_",k,".csv"))
  
  
  #### draw plot
  #CBCGAVStcga.gain.res = read.csv(paste0("/results/CBCGAVStcga_WHITE_gain_res_",k,".csv"),row.names = 1)
  p = ggplot(CBCGAVStcga.gain.res,aes(x=CBCGA, y = TCGA)) + 
    geom_point(alpha = 1, size = 2.5, col="#cecece") + 
    geom_point(data = CBCGAVStcga.gain.res[CBCGAVStcga.gain.res$adj.p < 0.05,] , 
               alpha = 1, size = 2.5, col="#C0052A") +
    geom_abline(intercept=0,slope=1, linetype="dashed", size = 1.2, color = '#747475')+
    xlim(0,1) + ylim(0,1) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA, colour = "black", size=1.2)) +
    geom_text_repel(data= CBCGAVStcga.gain.res[CBCGAVStcga.gain.res$adj.p < 0.05 & 
                                                 (CBCGAVStcga.gain.res$CBCGA > 0.25 | CBCGAVStcga.gain.res$TCGA > 0.25 ),] ,
                    aes(label=rownames(CBCGAVStcga.gain.res[CBCGAVStcga.gain.res$adj.p < 0.05  & 
                                                              (CBCGAVStcga.gain.res$CBCGA > 0.25 | CBCGAVStcga.gain.res$TCGA > 0.25 ),] )
                    ),
                    col="black",alpha = 1,size=3,
                    max.overlaps = 20) 
  ggsave(p,filename = paste0("/results/Fig3C_S3_CBCGAvsTCGA_WHITE_gain_",k,".pdf"),height = 4,width = 4)
  #export::graph2ppt(p,file = paste0("/results/CBCGAvsTCGA_WHITE_gain_",k,".pdf"),height = 4,width = 4)
  
  #CBCGAVStcga.loss.res = read.csv(paste0("/results/CBCGAVStcga_WHITE_loss_res_",k,".csv"),row.names = 1)
  p = ggplot(CBCGAVStcga.loss.res,aes(x=CBCGA, y = TCGA)) + 
    geom_point(alpha = 1, size = 2.5, col="#cecece") + 
    geom_point(data = CBCGAVStcga.loss.res[CBCGAVStcga.loss.res$adj.p < 0.05,] , 
               alpha = 1, size = 2.5, col="#4E85AC") +
    geom_abline(intercept=0,slope=1, linetype="dashed", size = 1.2, color = '#747475')+
    xlim(0,1) + ylim(0,1) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA, colour = "black", size=1.2)) +
    geom_text_repel(data= CBCGAVStcga.loss.res[CBCGAVStcga.loss.res$adj.p < 0.05 & 
                                                 (CBCGAVStcga.loss.res$CBCGA > 0.25 | CBCGAVStcga.loss.res$TCGA > 0.25 ),] ,
                    aes(label=rownames(CBCGAVStcga.loss.res[CBCGAVStcga.loss.res$adj.p < 0.05 &
                                                              (CBCGAVStcga.loss.res$CBCGA > 0.25 | CBCGAVStcga.loss.res$TCGA > 0.25 ),])),
                    col="black",alpha = 1,size=3,
                    max.overlaps = 20) 
  ggsave(p,filename = paste0("/results/Fig3C_S3_CBCGAvsTCGA_WHITE_loss_",k,".pdf"),height = 4,width = 4)
  #export::graph2ppt(p,file = paste0("/results/CBCGAvsTCGA_WHITE_loss_",k,".pdf"),height = 4,width = 4)
}


# Figure 2G-------------------------------------------------------------------
Tes_Cases <- CBCGA.Extended_Cohort.Info$PatientCode[CBCGA.Extended_Cohort.Info$CBCGA_Cohort == "Yes"]
HRpHER2p_Cases <- CBCGA.Extended_Cohort.Info$PatientCode[CBCGA.Extended_Cohort.Info$Clinical_Subtype == "HR+HER2+" ]

EXP_Tumor <- CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode[, str_detect(colnames(CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode), fixed("_T"))]
colnames(EXP_Tumor) <- substr( colnames(EXP_Tumor), 1, 4)
Tes_EXP   <- EXP_Tumor[, colnames(EXP_Tumor) %in% Tes_Cases]
Tes_EXP <- log2(Tes_EXP + 1)

Tes_PRO <- CBCGA.Extended_PRO_normalized
Tes_PRO <- Tes_PRO[, str_detect(colnames(Tes_PRO), fixed("_T"))]
colnames(Tes_PRO) <- substr( colnames(Tes_PRO), 1, 4)
Tes_PRO <- Tes_PRO[, colnames(Tes_PRO) %in% Tes_Cases]

Tes_CNA <- CBCGA_GISTICgene.alldata

### Generate Matrix
SCNA_Matrix <- Tes_CNA[c("NF1", "ACACA", "LASP1", "CDK12", "STARD3", "ERBB2", "GRB7", "MED24", "TOP2A"), ]
# SCNA_Matrix <- apply(SCNA_Matrix, 1, scale.default)
SCNA_Matrix[SCNA_Matrix > 5] <- 5

mRNA_Matrix <- Tes_EXP[c("NF1", "ACACA", "LASP1", "CDK12", "STARD3", "ERBB2", "GRB7",  "MED24", "TOP2A"), ]
mRNA_Matrix <- t(apply(t(mRNA_Matrix), 2, scale.default))
colnames(mRNA_Matrix) <- colnames(Tes_EXP)

PRO_Matrix  <- Tes_PRO[c("NF1", "ACACA", "LASP1", "CDK12", "STARD3", "ERBB2", "GRB7", "MED24", "TOP2A"), ]
PRO_Matrix  <- t(apply(t(PRO_Matrix), 2, scale.default))
colnames(PRO_Matrix) <- colnames(Tes_PRO)

Tes_Cases <- intersect(Tes_Cases, HRpHER2p_Cases)
Tes_Cases <- intersect( Tes_Cases, colnames(SCNA_Matrix) )
Tes_Cases <- intersect( Tes_Cases, colnames(mRNA_Matrix) )
Tes_Cases <- intersect( Tes_Cases, colnames(PRO_Matrix) ) 

Tes_CNA   <- SCNA_Matrix[, Tes_Cases]
Tes_Cases <- Tes_Cases[order(Tes_CNA["ERBB2", ], decreasing = T)]

Tes_CNA   <- Tes_CNA[, Tes_Cases]
Tes_EXP   <- EXP_Tumor[, Tes_Cases]
Tes_PRO   <- Tes_PRO[, Tes_Cases]

Clin_Matrix <- CBCGA.Extended_Cohort.Info[, c("PAM50_classifier")]

SCNA_Matrix <- SCNA_Matrix[, Tes_Cases]
mRNA_Matrix <- mRNA_Matrix[, Tes_Cases]
PRO_Matrix  <- PRO_Matrix[, Tes_Cases]

Pheatmap.breaks <- function(matrix, type, n.cut)
{
  pos    <- max(matrix) + 0.1 ; neg    <- min(matrix) + 0.1
  poscut <- n.cut
  negcut <- n.cut
  if(type == 1) 
  {mybreaks1 <- seq(neg, 0, (-neg)/negcut); mybreaks2 <- seq(0, pos, pos/poscut)[-1]
  } else {
    mybreaks1 <- -(seq(-sqrt(-neg), 0, sqrt(-neg)/negcut)) ^2
    mybreaks2 <- (seq(0, sqrt(pos), sqrt(pos)/poscut)[-1]) ^2
  } 
  mybreaks <- c(mybreaks1, mybreaks2)
  return(mybreaks)
}

Anno <- CBCGA_Cohort.Info[colnames(PRO_Matrix),  c("PR positive percentage",  "ER positive percentage", "Intrinsic subtype (PAM50)")]
colnames(Anno) <- c("PR", "ER", "PAM50")
Anno$PR[is.na(Anno$PR)] = 0
Anno$PR_percentage = Anno$PR
Anno$ER_percentage = Anno$ER

Anno$PR_percentage[Anno$PR == 0] = "0"
Anno$PR_percentage[Anno$PR < 10 & Anno$PR > 0 ] = "1-9"
Anno$PR_percentage[Anno$PR >= 10 & Anno$PR < 30] = "10-29"
Anno$PR_percentage[Anno$PR >= 30 & Anno$PR < 70] = "30-69"
Anno$PR_percentage[Anno$PR >= 70 ] = "70-100"

Anno$ER_percentage[Anno$ER < 10] = "1-9"
Anno$ER_percentage[Anno$ER >= 10 & Anno$ER < 30] = "10-29"
Anno$ER_percentage[Anno$ER >= 30 & Anno$ER < 70] = "30-69"
Anno$ER_percentage[Anno$ER >= 70 ] = "70-100"
Anno = Anno[,c("PR_percentage", "ER_percentage", "PAM50")]

Anno_Color <- list(
  PR_percentage =  c( "70-100" = "#000000", "30-69"  = "#666666","10-29" = "#999999" , "1-9" =  "#CCCCCC","0" = "#FFFFFF"), 
  ER_percentage = c( "70-100" = "#000000", "30-69"  = "#666666","10-29" = "#999999" , "1-9" =  "#CCCCCC"), 
  PAM50 = Color_PAM50
)


### Plots
## SCNA
n.cut        <- 50
mybreaks_CNA <- Pheatmap.breaks(SCNA_Matrix, 1, n.cut)
mycolor_CNA  <- c(colorRampPalette(c("#3A53A4", "#FFFFFF"))(n.cut), colorRampPalette(c( "#FFFFFF", "#E21E26"))(n.cut))

pheatmap(SCNA_Matrix, cluster_rows = F, cluster_cols = F, color = mycolor_CNA, breaks = mybreaks_CNA, border_color = NA, 
         annotation_col = Anno, annotation_colors = Anno_Color,
         filename = "/results/Fig3D_CNA.pdf", width = 10, height = 3.4)
## mRNA
n.cut        <- 50
mybreaks_mRNA <- Pheatmap.breaks(mRNA_Matrix, 2, n.cut)
mycolor_mRNA  <- c(colorRampPalette(c("Green", "#000000"))(n.cut), colorRampPalette(c( "#000000", "red"))(n.cut))

pheatmap(mRNA_Matrix, cluster_rows = F, cluster_cols = F, color = mycolor_mRNA, breaks = mybreaks_mRNA, border_color = NA,
         annotation_col = Anno, annotation_colors = Anno_Color,
         filename = "/results/Fig3D_mRNA.pdf", width = 10, height = 3.4)


## Protein
n.cut        <- 50
mybreaks_Pro <- Pheatmap.breaks(PRO_Matrix, 1, n.cut)
mycolor_Pro  <- c(colorRampPalette(c("#3A53A4", "#FFFFFF"))(n.cut), colorRampPalette(c( "#FFFFFF", "#E21E26"))(n.cut))
pheatmap(PRO_Matrix, cluster_rows = F, cluster_cols = F, color = mycolor_Pro, breaks = mybreaks_Pro, border_color = NA,
         annotation_col = Anno, annotation_colors = Anno_Color,
         filename = "/results/Fig3D_Pro.pdf", width = 10, height = 3.4)