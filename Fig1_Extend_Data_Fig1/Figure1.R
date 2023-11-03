rm(list = ls()) ; graphics.off() ; gc()
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(export)
library(ggsci)
library(RColorBrewer)
library(dendextend)
load("Fig1_Supp_data.Rdata")
load("col_order230424.RData")
load("CBCGA.Extended_MergedData_V2.6_220802.Rdata")

################ 1. prepare clinical data ---------
CBCGAClin = CBCGA.Extended_Cohort.Info[CBCGA.Extended_Cohort.Info$CBCGA_Cohort == "Yes",]
CBCGAClin$TMB = TMB[rownames(CBCGAClin)]

CBCGAClin$Lymph_node_status = LN[rownames(CBCGAClin)]
CBCGAClin$Lymph_node_status[CBCGAClin$Lymph_node_status > 0] = 1

CBCGAClin$Ki67 = Ki67[rownames(CBCGAClin)]

CBCGAClin = CBCGAClin[colnames(CBCGA_merge_mat),]

heatmap_order <- col_order

CBCGAClin <- CBCGAClin[match(col_order, CBCGAClin$PatientCode), ]
CBCGAClin$Pathomics = CBCGA_Cohort.Info[rownames(CBCGAClin),"Pathomics"]
CBCGAClin$Radiomics = CBCGA_Cohort.Info[rownames(CBCGAClin),"Radiomics"]
CBCGA_merge_mat <- CBCGA_merge_mat[, col_order]

################ 2. prepare polar and lipid data ---------

load("res_diff_metaPol_vsOther.Rdata")

# polar metabolite
DEG_sig <- bind_rows(res_diff)
DEG_sig <- DEG_sig[DEG_sig$Sig == "Up", ]
DEG_sig <- group_by(DEG_sig, Sig)%>%
  arrange(desc(abs(logFC)))%>%
  filter(abs(logFC) > 1.5)
Picked_polar <- unique(DEG_sig$Symbol)

# lipids
load("res_diff_metaLipid_vsOther_PAM50.Rdata")
DEG_sig <- bind_rows(res_diff)
DEG_sig <- DEG_sig[DEG_sig$Sig == "Up", ]
DEG_sig <- group_by(DEG_sig, Sig)%>%
  arrange(desc(abs(logFC)))%>%
  filter(abs(logFC) > 1.5)
Picked_lipid <- unique(DEG_sig$Symbol)


data <- CBCGA.Extended_pol
info <- CBCGA_pol_anno

id_dat <- names(data)
id_dat <- id_dat[str_detect(id_dat, "_T$")]
id_dat <- id_dat[substring(id_dat, 1, 4)%in%sample_ord]


info_polar <- info[info$peak%in%Picked_polar, ];table(info_polar$Metabolite_class)
dat_polar <- data[, id_dat];names(dat_polar) <- substring(names(dat_polar), 1, 4)


tmp = as.data.frame( matrix(nrow = nrow(dat_polar), 
                            ncol = length(setdiff(sample_ord, colnames(dat_polar)) ) ) )
rownames(tmp) = rownames(dat_polar); colnames(tmp) = setdiff(sample_ord, colnames(dat_polar))

dat_polar_addNA <- cbind(dat_polar, tmp)

# 
CBCGAPolar_mat <- dat_polar_addNA[Picked_polar, heatmap_order]
identical(heatmap_order, names(CBCGAPolar_mat))


CBCGAPolar_mat_colname <- colnames(CBCGAPolar_mat)
CBCGAPolar_mat <- apply(CBCGAPolar_mat, 1, scale) %>% t()
colnames(CBCGAPolar_mat) <- CBCGAPolar_mat_colname

range(CBCGAPolar_mat, na.rm = T)

#
data <- CBCGA.Extended_lip
info <- CBCGA_lip_anno


id_dat <- names(data)
id_dat <- id_dat[str_detect(id_dat, "_T$")]
id_dat <- id_dat[substring(id_dat, 1, 4)%in%sample_ord]

# 
info_lipid <- info[rownames(info)%in%Picked_lipid, ];table(info_lipid$Lipid.super.class)
dat_lipid <- data[, id_dat];names(dat_lipid) <- substring(names(dat_lipid), 1, 4)

# 
tmp = as.data.frame( matrix(nrow = nrow(dat_lipid), 
                            ncol = length(setdiff(sample_ord, colnames(dat_lipid)) ) ) )
rownames(tmp) = rownames(dat_lipid); colnames(tmp) = setdiff(sample_ord, colnames(dat_lipid))

dat_lipid_addNA <- cbind(dat_lipid, tmp)

# 
CBCGAlipid_mat <- dat_lipid_addNA[Picked_lipid, heatmap_order]
identical(heatmap_order, names(CBCGAlipid_mat))#T

CBCGAlipid_mat_colname <- colnames(CBCGAlipid_mat)
CBCGAlipid_mat <- apply(CBCGAlipid_mat, 1, scale) %>% t()
colnames(CBCGAlipid_mat) <- CBCGAlipid_mat_colname

range(CBCGAlipid_mat, na.rm = T)

################ 3. prepare RNA data ---------
if(TRUE){
  data <- CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode
  
  id_dat <- names(data)
  id_dat <- id_dat[str_detect(id_dat, "_T$")]
  id_dat <- id_dat[substring(id_dat, 1, 4)%in%sample_ord]
  
  dat_RNA<- data[, id_dat];names(dat_RNA) <- substring(names(dat_RNA), 1, 4)
  
  
  tmp <- as.data.frame(matrix(nrow = nrow(dat_RNA), 
                              ncol = length(setdiff(sample_ord, colnames(dat_RNA)) ) ) )
  rownames(tmp) <- rownames(dat_RNA) ; colnames(tmp) = setdiff(sample_ord, colnames(dat_RNA))
  
  dat_RNA_addNA <- log2( cbind(dat_RNA,tmp) +1 )
  
  # 
  CBCGARNA_mat <- dat_RNA_addNA[rownames(picked.RNA), heatmap_order]
  identical(heatmap_order, names(CBCGARNA_mat))
  
  
  CBCGARNA_mat_colname <- colnames(CBCGARNA_mat)
  CBCGARNA_mat <- apply(CBCGARNA_mat, 1, scale) %>% t()
  colnames(CBCGARNA_mat) <- CBCGARNA_mat_colname
  
  range(CBCGARNA_mat, na.rm = T)

}

################ 4. prepare Protein data ---------
if(TRUE){
  data <- CBCGA.Extended_PRO_normalized
  
  id_dat <- names(data)
  id_dat <- id_dat[str_detect(id_dat, "_T$")]
  id_dat <- id_dat[substring(id_dat, 1, 4)%in%sample_ord]
  
  dat_PRO <- data[, id_dat];names(dat_PRO) <- substring(names(dat_PRO), 1, 4)
  
  tmp <- as.data.frame(matrix(nrow = nrow(dat_PRO), 
                              ncol = length(setdiff(sample_ord, colnames(dat_PRO)) ) ) )
  rownames(tmp) <- rownames(dat_PRO) ; colnames(tmp) = setdiff(sample_ord, colnames(dat_PRO))
  
  dat_PRO_addNA <- cbind(dat_PRO,tmp)
  
  # 
  CBCGAPRO_mat <- dat_PRO_addNA[rownames(picked.Pro), heatmap_order]
  identical(heatmap_order, names(CBCGAPRO_mat))#T
  
  CBCGAPRO_mat_colname <- colnames(CBCGAPRO_mat)
  CBCGAPRO_mat <- apply(CBCGAPRO_mat, 1, scale) %>% t()
  colnames(CBCGAPRO_mat) <- CBCGAPRO_mat_colname
  
  range(CBCGAPRO_mat, na.rm = T)#-6 to 7

}

################ 5. prepare annotation ---------
all(rownames(picked.Pro) == rownames(CBCGAPRO_mat))
all(rownames(picked.RNA) == rownames(CBCGARNA_mat))

# RNA
picked.RNA$path <- factor(picked.RNA$path, levels = c('ESTROGEN_RESPONSE_EARLY', 'ESTROGEN_RESPONSE_LATE', 'ERBB2_SIGNALING_PATHWAY',
                                                      'SIGNALING_BY_ERBB2_ECD_MUTANTS', 'NATURAL_KILLER_CELL_MEDIATED_IMMUNITY',
                                                      'POSITIVE_REGULATION_OF_INNATE_IMMUNE_RESPONSE', 'T_CELL_RECEPTOR_SIGNALING_PATHWAY',
                                                      'MYC_TARGETS', 'CELL_CYCLE', 'P53_SIGNALING_PATHWAY', 'NEGATIVE_REGULATION_OF_DNA_REPAIR'))
picked.RNA$name <- rownames(picked.RNA)
picked.RNA <- picked.RNA[order(picked.RNA$path), ]
picked.RNA <- picked.RNA[, -2, drop = F]

# Protein
picked.Pro$path <- factor(picked.Pro$path, levels = c('ESTROGEN_RESPONSE_EARLY', 'ESTROGEN_RESPONSE_LATE', 'ERBB2_SIGNALING_PATHWAY',
                                                      'SIGNALING_BY_ERBB2_ECD_MUTANTS', 'NATURAL_KILLER_CELL_MEDIATED_IMMUNITY',
                                                      'POSITIVE_REGULATION_OF_INNATE_IMMUNE_RESPONSE', 'T_CELL_RECEPTOR_SIGNALING_PATHWAY',
                                                      'MYC_TARGETS', 'CELL_CYCLE', 'P53_SIGNALING_PATHWAY', 'NEGATIVE_REGULATION_OF_DNA_REPAIR'))
picked.Pro$name <- rownames(picked.Pro)
picked.Pro <- picked.Pro[order(picked.Pro$path), ]
picked.Pro <- picked.Pro[, -2, drop = F]

right_annotation_RNA <- rowAnnotation(df = picked.RNA,
                                      show_annotation_name = F,
                                      col = list(path = c("T_CELL_RECEPTOR_SIGNALING_PATHWAY" = brewer.pal(12, "Set3")[1],
                                                          "POSITIVE_REGULATION_OF_INNATE_IMMUNE_RESPONSE" = brewer.pal(12, "Set3")[2] ,
                                                          "NATURAL_KILLER_CELL_MEDIATED_IMMUNITY" = brewer.pal(12, "Set3")[3],
                                                          "INNATE_IMMUNE_RESPONSE_ACTIVATING_SIGNAL_TRANSDUCTION" = brewer.pal(12, "Set3")[4],
                                                          "SIGNALING_BY_ERBB2_ECD_MUTANTS" = brewer.pal(12, "Set3")[5],
                                                          "ERBB2_SIGNALING_PATHWAY" = brewer.pal(12, "Set3")[6],
                                                          "ESTROGEN_RESPONSE_EARLY" = brewer.pal(12, "Set3")[7],
                                                          "ESTROGEN_RESPONSE_LATE" = brewer.pal(12, "Set3")[8],
                                                          "MYC_TARGETS" = brewer.pal(12, "Set3")[9],
                                                          "CELL_CYCLE" = brewer.pal(12, "Set3")[10],
                                                          "P53_SIGNALING_PATHWAY" = brewer.pal(12, "Set3")[11],
                                                          "NEGATIVE_REGULATION_OF_DNA_REPAIR" = brewer.pal(12, "Set3")[12]))
)

right_annotation_Pro <- rowAnnotation(df = picked.Pro,
                                      show_annotation_name = F,
                                      col = list(path = c("T_CELL_RECEPTOR_SIGNALING_PATHWAY" = brewer.pal(12, "Set3")[1],
                                                          "POSITIVE_REGULATION_OF_INNATE_IMMUNE_RESPONSE" = brewer.pal(12, "Set3")[2] ,
                                                          "NATURAL_KILLER_CELL_MEDIATED_IMMUNITY" = brewer.pal(12, "Set3")[3],
                                                          "INNATE_IMMUNE_RESPONSE_ACTIVATING_SIGNAL_TRANSDUCTION" = brewer.pal(12, "Set3")[4],
                                                          "SIGNALING_BY_ERBB2_ECD_MUTANTS" = brewer.pal(12, "Set3")[5],
                                                          "ERBB2_SIGNALING_PATHWAY" = brewer.pal(12, "Set3")[6],
                                                          "ESTROGEN_RESPONSE_EARLY" = brewer.pal(12, "Set3")[7],
                                                          "ESTROGEN_RESPONSE_LATE" = brewer.pal(12, "Set3")[8],
                                                          "MYC_TARGETS" = brewer.pal(12, "Set3")[9],
                                                          "CELL_CYCLE" = brewer.pal(12, "Set3")[10],
                                                          "P53_SIGNALING_PATHWAY" = brewer.pal(12, "Set3")[11],
                                                          "NEGATIVE_REGULATION_OF_DNA_REPAIR" = brewer.pal(12, "Set3")[12]))
                                      
)

CBCGARNA_mat <- CBCGARNA_mat[rownames(picked.RNA), ]
CBCGAPRO_mat <- CBCGAPRO_mat[rownames(picked.Pro), ]

tmp = data.frame(row.names = c(rownames(CBCGAPolar_mat)),
                 class = c(info_polar$Metabolite_class[match(rownames(CBCGAPolar_mat), info_polar$peak)]))
color_polar <-  structure(c("#8FC3E0", "#95C88A", "#E51C2C", "#593D8E", "#0A74B3", "#F9B466", "#A6A9A9", "#B9A4C4"), 
                          names = c("Amino acid", "Carbohydrates", "Lipid", "Nucleotide", "Peptide", "Vitamins and Cofactors",
                                    "Xenobiotics", "Other"))

tmp_polar <- tmp

names(tmp) <- "Polar_Class"
right_annotation_polar = rowAnnotation(df = tmp,
                                       show_annotation_name = F,
                                       col = list(Polar_Class = color_polar[unique(tmp$Polar_Class)])
)

tmp = data.frame(row.names = c(rownames(CBCGAlipid_mat)),
                 class = c(info_lipid$Lipid.super.class[match(rownames(CBCGAlipid_mat), rownames(info_lipid))]))
color_lipid <- structure(c("#E82029", "#1F78B6", "#36A754", "#593D8E"), 
                         names = c("Fatty acyls [FA]",  "Glycerolipids [GL]", "Glycerophospholipids [GP]", "Sphingolipids [SP]"))

tmp_lipid <- tmp

names(tmp) <- "Lipid_Class"
right_annotation_lipid = rowAnnotation(df = tmp,
                                       show_annotation_name = F,
                                       col = list(Lipid_Class = color_lipid)
)


################ 6. prepare color ---------

if(TRUE){
  col = c("Splicing" = "#6b4694", "Missense" = "#4d93cd", "Nonsense" = "#e63a3a", "Frameshift" = '#fab71b', "Inframe" = "#ecd71e",
          "Nonstop" = '#ab5b9e', "Translation_Start_Site" = '#018b38', "No" = '#eaeaea', 'Amplification' = '#BC102B', 'Gain' = '#DF989E',
          'Loss' = '#AEC4D6', 'Deletion' = '#5385AC')
  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.5, "pt"), h-unit(2.5, "pt"),
                gp = gpar(fill = "#CCCCCC", col = NA))
    },
    # small purple
    Splicing = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.5, "pt"), h*0.67,
                gp = gpar(fill = col["Splicing"], col = NA))
    },
    # small green
    Missense = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.5, "pt"), h*0.33,
                gp = gpar(fill = col["Missense"], col = NA))
    },
    # small black
    Nonsense = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.5, "pt"), h*0.33,
                gp = gpar(fill = col["Nonsense"], col = NA))
    },
    # small black
    Frameshift = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.5, "pt"), h*0.33,
                gp = gpar(fill = col["Frameshift"], col = NA))
    },
    # small black
    Inframe = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.5, "pt"), h*0.33,
                gp = gpar(fill = col["Inframe"], col = NA))
    },
    # small black
    Nonstop = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.5, "pt"), h*0.33,
                gp = gpar(fill = col["Nonstop"], col = NA))
    },
    # small black
    Translation_Start_Site = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.5, "pt"), h*0.33,
                gp = gpar(fill = col["Translation_Start_Site"], col = NA))
    },
    Amplification = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.5, "pt"), h-unit(2.5, "pt"),
                gp = gpar(fill = col["Amplification"], col = NA))
    },
    Gain = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.5, "pt"), h-unit(2.5, "pt"),
                gp = gpar(fill = col["Gain"], col = NA))
    },
    Loss = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.5, "pt"), h-unit(2.5, "pt"),
                gp = gpar(fill = col["Loss"], col = NA))
    },
    Deletion = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.5, "pt"), h-unit(2.5, "pt"),
                gp = gpar(fill = col["Deletion"], col = NA))
    },
    No = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.5, "pt"), h-unit(2.5, "pt"),
                gp = gpar(fill = col["No"], col = NA))
    }
  )
  
  column_title <- ""
  heatmap_legend_param <- list(title = "", at = c("Splicing", "Missense", 'Nonsense', 'Frameshift', 'Inframe', "Nonstop", "Translation_Start_Site",
                                                  'Amplification', 'Gain', 'Loss', 'Deletion', 'No'),
                               labels = c("Splicing", "Missense", 'Nonsense', 'Frameshift', 'Inframe', "Nonstop", "Translation_Start_Site",
                                          'Amplification', 'Gain', 'Loss', 'Deletion', 'No'))
}

## MutSig
mutsig_ord = c("SBS5","SBS1",
               "SBS13","SBS2",
               "SBS8","SBS18","SBS17a","SBS17b",
               "SBS20","SBS26","SBS30","SBS3","SBS6",
               "SBSNA")
mutsig_col = c('SBS5' ='#0084C8' , 'SBS1' = '#A6CEE3', 
               'SBS13' =  '#009100' ,'SBS2' = '#9ADE00' ,
               'SBS8' ='#6A3D9A'  ,
               'SBS17a' ='#33C6BB', 'SBS17b' =  '#99E3DD', 'SBS18' = '#008A80',
               'SBS3' = '#FF6600','SBS6' = '#DC0000', 'SBS20' = '#FFFF3E','SBS26' = '#FFC022','SBS30' = '#FF9900', 
               'SBSNA' = '#eaeaea')

MutSig = MutSig[mutsig_ord,]
mutsig_col = mutsig_col[mutsig_ord]
MutSig <- MutSig[, match(heatmap_order, colnames(MutSig))]

col_fun_age <- colorRamp2(c(24, 53, 90), c('#ecf3e4', '#9ec27d', '#489205'))
col_fun_ki67 <- colorRamp2(c(3, 30, 98), c('#EBEFF6', '#99AFD2', '#3C74AE'))
col_fun_hrd <- colorRamp2(c(0, 16, 75), c('#FDF9E7', '#EFE089', '#D7C801'))

pos <- 5 + 0.1
neg <- -5 + 0.1
poscut <- 50
negcut <- 50
mybreaks1 <- -(seq(-sqrt(-neg), 0, sqrt(-neg)/negcut)) ^2
mybreaks2 <- (seq(0, sqrt(pos), sqrt(pos)/poscut)[-1]) ^2
mybreaks  <- c(mybreaks1, mybreaks2)

# TMB
TMB1 <- as.numeric(CBCGAClin$TMB)
TMB1[which(TMB1 > 6)] <- NA
TMB2 <- as.numeric(CBCGAClin$TMB)
TMB2[which(TMB2 <= 6)] <- NA

################ 7. plot ---------
h_main <- oncoPrint(CBCGA_merge_mat, 
                    alter_fun = alter_fun, 
                    col = col, 
                    show_pct = T, 
                    pct_side = 'right',
                    column_title = '', 
                    heatmap_legend_param = heatmap_legend_param, 
                    column_order = heatmap_order,
                    row_names_side = "left", 
                    show_column_names = F, 
                    remove_empty_columns = F,
                    row_names_gp = gpar(fontsize = 9), 
                    column_title_gp = gpar(fontsize = 0),
                    column_split = factor(CBCGAClin$PAM50_classifier,
                                          levels = c("LumA", "LumB", "Her2", "Basal", "Normal")),
                    row_split = factor(c(rep('Somatic', 13),
                                         rep('Germline', 5),
                                         rep('AMP', 9), 
                                         rep('DEL', 8)), 
                                       levels = c('Somatic', 'Germline', 'AMP', 'DEL')),
                    top_annotation = HeatmapAnnotation(TMB2 = anno_points(TMB2, 
                                                                          size = unit(1, 'mm'),
                                                                          border = T,
                                                                          ylim = c(6, 20),
                                                                          height = unit(0.75, 'cm'),
                                                                          axis_param = list(at = seq(8, 20, by = 6),
                                                                                            labels = seq(8, 20, by = 6))),
                                                       TMB1 = anno_points(TMB1,
                                                                          size = unit(1, 'mm'), 
                                                                          border = T, 
                                                                          ylim = c(0, 5.95), 
                                                                          height = unit(1.5, 'cm'),
                                                                          axis_param = list(at = c(0, 2, 4, 6),
                                                                                            labels = c(0, 2, 4, 6))),
                                                       Intrinsic_subtype = CBCGAClin$PAM50_classifier,
                                                       IHC_Subtype = CBCGAClin$Clinical_Subtype,
                                                       Age.at.surgery = CBCGAClin$Age,
                                                       Lymph_node_status = CBCGAClin$Lymph_node_status,
                                                       Ki67 = as.numeric(CBCGAClin$Ki67),
                                                       HRD = as.numeric(CBCGAClin$HRD),
                                                       Relapse = factor(as.character(CBCGAClin$RFS_status)),
                                                       Pathomics = CBCGAClin$Pathomics,
                                                       Radiomics = CBCGAClin$Radiomics,
                                                       COSMIC_sig = anno_barplot(t(MutSig), height = unit(1.5, 'cm'),
                                                                                 gp = gpar(fill = mutsig_col,  col = NA)),
                                                       annotation_name_side = "left",
                                                       annotation_name_gp = gpar(fontsize = 9),
                                                       col = list(IHC_Subtype = c('HR+HER2-' = '#0085c4', 'HR+HER2+' = '#7ab801',
                                                                                  'HR-HER2+' = '#f2af01', 'TNBC' = '#dc5035'),
                                                                  Intrinsic_subtype = c('LumA' = '#0077c9', 'LumB' = '#74d2e8',
                                                                                        'Her2' = '#7552cd', 'Basal' = '#e4002c',
                                                                                        'Normal' = '#cecece',"Unknown" = "grey"),
                                                                  Age.at.surgery = col_fun_age,
                                                                  Lymph_node_status = c('0' = 'white', '1' = 'black'),
                                                                  Ki67 = col_fun_ki67,
                                                                  HRD = col_fun_hrd,
                                                                  Relapse = c('0' = 'white', '1' = 'black','NA' = "grey"),
                                                                  Pathomics = c('No' = 'white', 'Yes' = '#6A4984'),
                                                                  Radiomics = c('No' = 'white', 'Yes' = '#D27683'),
                                                                  COSMIC_sig = mutsig_col)),
                    gap = unit(0.1, 'inch'))

# RNA
h_RNA <- Heatmap(CBCGARNA_mat, 
                 name = "Median centered \nlog2(FPKM+1)",
                 show_row_names = F,
                 show_column_names = F, 
                 column_order = heatmap_order,
                 column_split = factor(CBCGAClin$PAM50_classifier, 
                                       levels = c("LumA", "LumB", "Her2", "Basal", "Normal")),
                 column_title = NULL,
                 height = unit(4, 'cm'),
                 cluster_columns = F, 
                 cluster_rows = F, 
                 na_col = '#eaeaea',
                 col = colorRamp2(mybreaks, c(colorRampPalette(c("green", "#000000"))(51), 
                                              colorRampPalette(c("#000000", "red"))(50))), 
                 right_annotation = right_annotation_RNA) 


h_protein <- Heatmap(CBCGAPRO_mat, 
                     name = "Normlized \nabundance of proteins",
                     show_row_names = F,
                     show_column_names = F, 
                     column_order = heatmap_order,
                     column_split = factor(CBCGAClin$PAM50_classifier, 
                                           levels = c("LumA", "LumB", "Her2", "Basal", "Normal")),
                     column_title = NULL,
                     height = unit(4, 'cm'),
                     cluster_columns = F,
                     cluster_rows = F, 
                     na_col = '#eaeaea',
                     col = colorRamp2(mybreaks, c(colorRampPalette(c("#3A53A4", "#FFFFFF"))(51),
                                                  colorRampPalette(c( "#FFFFFF", "#E21E26"))(50))), 
                     right_annotation = right_annotation_Pro)


# polar
gg <- hclust(dist(CBCGAPolar_mat))
dend <- as.dendrogram(gg)
row_clust <- dend%>%ladderize(FALSE)%>%as.hclust()

h_polar <- Heatmap(CBCGAPolar_mat, 
                   name = "Normlized \nabundance of metabolites",
                   show_row_names = F, 
                   # row_names_gp = gpar(fontsize = 5),
                   show_column_names = F, 
                   height = unit(2, 'cm'),
                   cluster_columns = F, 
                   cluster_rows = row_clust,
                   # row_split = tmp_polar$class,
                   column_order = heatmap_order,
                   column_split = factor(CBCGAClin$PAM50_classifier, 
                                         levels = c("LumA", "LumB", "Her2", "Basal", "Normal")),
                   column_title = NULL,
                   row_gap = unit(0, "cm"),
                   row_title = NULL,
                   na_col = '#eaeaea',
                   show_row_dend = F,
                   # top_annotation = topanno,
                   right_annotation = right_annotation_polar,
                   # col = colorRamp2(c(-5, 0, 5), c("#26B2E3", "#000000", "#F4EB08"))
                   col = colorRamp2(mybreaks, c(colorRampPalette(c("#26B2E3", "#000000"))(51), 
                                                colorRampPalette(c( "#000000", "#F4EB08"))(50)))
)


# lipid
gg <- hclust(dist(CBCGAlipid_mat))
dend <- as.dendrogram(gg)
# dend <- reorder(dend, wts = )
row_clust_lipid <- dend%>%ladderize(FALSE)%>%as.hclust()

h_lipid <- Heatmap(CBCGAlipid_mat, 
                   name = "Normlized \nabundance of metabolites",
                   show_row_names = F, 
                   show_column_names = F, 
                   # row_names_gp = gpar(fontsize = 5),
                   height = unit(2, 'cm'),
                   cluster_columns = F, 
                   cluster_rows = row_clust_lipid,
                   # row_split = tmp_lipid$class,
                   column_order = heatmap_order,
                   column_split = factor(CBCGAClin$PAM50_classifier, 
                                         levels = c("LumA", "LumB", "Her2", "Basal", "Normal")),
                   column_title = NULL,
                   row_gap = unit(0, "cm"),
                   row_title = NULL,
                   na_col = '#eaeaea',
                   show_row_dend = F,
                   right_annotation = right_annotation_lipid,
                   # col = colorRamp2(c(-5, 0, 5), c("#26B2E3", "#000000", "#F4EB08"))
                   col = colorRamp2(mybreaks, c(colorRampPalette(c("#26B2E3", "#000000"))(51), 
                                                colorRampPalette(c( "#000000", "#F4EB08"))(50)))
)


CBCGA_landscape <- h_main%v%h_RNA%v%h_protein%v%h_polar%v%h_lipid

# plot
pdf(file = "/Volumes/零号机/研究生科研/科研项目/已完成的项目/研一/CBCGA/接收后修改/Fig1_FigS1_TableS1自查/Figure1_CBCGA_landscape_pam50.pdf",
    width = 18, height = 18)
CBCGA_landscape

dev.off()




