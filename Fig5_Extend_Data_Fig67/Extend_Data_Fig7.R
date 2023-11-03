###########################################  Fig.S7 Fusion genes in CBCGA ############################################
rm(list = ls())
#### ----------------library---------------- ####
library(tidyverse)
#library(ggpubr)
library(export)
library(foreach)
library(reshape2)
library(circlize)
#library(readxl)
library(survival)
#library(survminer)

#### ----------------import data---------------- ####
load('/data/FigS7_data.Rdata')

#### ----------------data prepare---------------- ####
TT_id <- meta_fusion_CBCGA$Gene_fusion_ID[meta_fusion_CBCGA$group=='T']
TP_id <- meta_fusion_CBCGA$Gene_fusion_ID[meta_fusion_CBCGA$group=='N']
meta_TT <- meta_fusion_CBCGA[meta_fusion_CBCGA$group=='T',]
meta_TP <- meta_fusion_CBCGA[meta_fusion_CBCGA$group=='N',]

#### ----------------Fig.S7 A---------------- ####
# gf_TT <- Fusion_gene_pair_ProteinCoding[,TT_id]
# gf_TT$Gene_pair <- rownames(gf_TT)
# gf_TT <- left_join(gf_TT,Fusion_transcript_ProteinCoding_detailed_info[,c('Gene_pair','Gene1','Gene1_chr','Gene2','Gene2_chr')],by='Gene_pair')
# 
# ch <- unique(c(gf_TT$Gene1_chr,gf_TT$Gene2_chr))
# ch <- factor(ch,levels = c("chr1", "chr2", "chr3" , "chr4" , "chr5"  ,"chr6" , "chr7" , "chr8" , "chr9" , "chr10" ,"chr11", "chr12", "chr13" ,"chr14" ,"chr15" ,"chr16", "chr17" ,"chr18" ,"chr19" , "chr20", "chr21", "chr22","chrX",'chrY'))
# ch <- factor(ch,levels = rev(levels(ch)))
# 
# gf_TT_ch <- data.frame()
# for (i in 1:24) {
#   df <- gf_TT[gf_TT$Gene1_chr==levels(ch)[i]|gf_TT$Gene2_chr==levels(ch)[i],]
#   gf_TT_ch <- rbind(gf_TT_ch,apply(df[,c(1:(ncol(gf_TT)-5))], 2, function(x){sum(x!=0)}))
# }
# gf_TT_ch[gf_TT_ch!=0] <- 1
# colnames(gf_TT_ch) <- colnames(gf_TT[,1:(ncol(gf_TT)-5)])
# rownames(gf_TT_ch) <- levels(ch)
# gf_TT_ch$ch <- rownames(gf_TT_ch)
# gf_TT_ch_reshape <- melt(gf_TT_ch,varing=list(1:(ncol(gf_TT)-5)),value.name = 'fusion')
# colnames(gf_TT_ch_reshape)[2] <- 'Gene_fusion_ID'
# gf_TT_ch_reshape <- left_join(gf_TT_ch_reshape,meta_TT[,c('Gene_fusion_ID','Clinical_Subtype','PAM50_classifier')],by='Gene_fusion_ID')
# gf_TT_ch_reshape$Clinical_Subtype <- factor(gf_TT_ch_reshape$Clinical_Subtype,levels =c('HR+HER2-','HR+HER2+','HR-HER2+','TNBC'))
# 
# 
# df_melt <- gf_TT_ch_reshape
# cols <- c('HR+HER2-'='#0085c4','HR+HER2+'='#7ab801','HR-HER2+'='#f2af01','TNBC'='#dc5035')                  
# p <- ggbarplot(df_melt,x="ch",y="fusion",fill="Clinical_Subtype",color = NA,
#                stat = 'identity')+
#   theme_classic()+
#   theme(panel.grid=element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(colour = 'black',size = 10),
#         axis.text.y = element_text(colour = 'black',size = 10),
#         legend.title = element_text(colour = 'black',size = 10))+
#   labs(x='',y='Samples with gene fusion event',title = 'Distribution of fusion genes across chromosomes',fill='Clinical subtype')+
#   scale_fill_manual(values = cols)+coord_flip()
# p
# #export::graph2ppt(p,'Fig_S7A_Distribution of fusion genes across chromosomes.ppt',height = 5, width = 5,append = F)
# ggsave(p,filename = '/results/Fig_S7A_Distribution of fusion genes across chromosomes.pdf',p,width = 5,height = 5)



#### ----------------Fig.S7 B---------------- ####
## Fusion genes in TP
gf_TP <- Fusion_gene_pair_ProteinCoding[,TP_id]
exist_in_sample <- apply(gf_TP, 1, function(x){sum(x!=0)})
FusionGene_TP <- names(exist_in_sample[exist_in_sample!=0])

## Fusion genes in TT: at least in more than two patients
gf_TT <- Fusion_gene_pair_ProteinCoding[,TT_id]
exist_in_sample <- apply(gf_TT, 1, function(x){sum(x!=0)})
exist_in_sample <- sort(exist_in_sample,decreasing = T)
FusionGene_TT <- setdiff(names(exist_in_sample[exist_in_sample>2]),FusionGene_TP)

## parameters
Fusion <- gf_TT[FusionGene_TT,]
Fusion$exist_in_sample <- apply(Fusion, 1, function(x){sum(x!=0)})
Fusion$Gene_pair <- rownames(Fusion)
Fusion$Gene1 <- foreach(i=1:nrow(Fusion),.combine=c) %do% strsplit(Fusion$Gene_pair[i],split="--")[[1]][1]
Fusion$Gene2 <- foreach(i=1:nrow(Fusion),.combine=c) %do% strsplit(Fusion$Gene_pair[i],split="--")[[1]][2]
Fusion <- Fusion[,c('Gene1','Gene2','Gene_pair','exist_in_sample')]
fs <- unique(c(Fusion$Gene1,Fusion$Gene2))
fs_len=foreach(i=1:length(fs),.combine = c) %do% max(Fusion$exist_in_sample[Fusion$Gene1==fs[i]|Fusion$Gene2==fs[i]]) 
names(fs_len)=fs 
fs_col=rep("red",length(fs)) 
names(fs_col)=fs
#database <- read_xlsx('Fusion_genes_database.xlsx')
fg_reported <- database$Gene_pair[database$database_merge==1]
fs_col[unique(unlist(strsplit(fg_reported,split="--")))]="black"  
factors=rep(fs,times=fs_len)

## draw 
pdf('/results/FigS7B_Recurrent_Fusion_filter_TP_rankByLetter.pdf',height=10,width=10)
set.seed(13111)
par(mfrow=c(1,1),mar=rep(8,4))
circos.par(cell.padding = c(0, 0, 0, 0),
           start.degree =90,  
           gap.degree =0.5, 
           track.margin = c(0, 0.02)) 
circos.initialize(factors, 
                  xlim = cbind(rep(0,length(fs)),fs_len)) 
circos.track(ylim = c(0, 1), 
             track.height = 0.05, 
             bg.col = rand_color(length(fs),luminosity = "random"), 
             bg.border = NA,
             panel.fun=function(x,y){
               sector.index = CELL_META$sector.index
               m = fs_len[sector.index]
               circos.text(m/2,1.5, 
                           label=names(m), 
                           col=fs_col[sector.index], 
                           cex=0.8, 
                           niceFacing=T, 
                           facing="reverse.clockwise",
                           adj=c(1,0.5),
                           xpd=T)
             }) 
for(i in 1:nrow(Fusion)) {
  circos.link(as.character(Fusion[i,1]),
              c(0,Fusion[i,4]), 
              as.character(Fusion[i,2]), 
              c(0,Fusion[i,4]), 
              col = rand_color(1, transparency = 0.4), border = NA)
}
title('Recurrent_Fusion_filtered_TP(n>2)')
circos.clear()
graphics.off()


#### ----------------Fig.S7 C---------------- ####
# df <- gf_TT[FusionGene_TT,] %>% t() %>% as.data.frame()
# df[df!=0] <- 1
# df$Gene_fusion_ID <- rownames(df)
# df <- merge(meta_TT,df,by='Gene_fusion_ID')
# df <- df[,c('Gene_fusion_ID','Clinical_Subtype',FusionGene_TT)]
# 
# df_melt <- melt(df,varying = list(fusionGene_TT),variable.name = 'Gene_pair',value.name = 'Fusion')
# df_melt$Gene_pair <- fct_inorder(df_melt$Gene_pair)
# df_melt$Gene_pair <- factor(df_melt$Gene_pair,levels = rev(levels(df_melt$Gene_pair)))
# df_melt$Clinical_Subtype <- fct_relevel(df_melt$Clinical_Subtype,c('HR+HER2-','HR+HER2+','HR-HER2+','TNBC'))
# 
# cols <- c('HR+HER2-'='#0085c4','HR+HER2+'='#7ab801','HR-HER2+'='#f2af01','TNBC'='#dc5035')
# p <- ggbarplot(df_melt,x="Gene_pair",y="Fusion",fill="Clinical_Subtype",color = NA,
#                stat = 'identity')+
#   border(color=NA,size = 0.1)+
#   theme_classic()+
#   theme(panel.grid=element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(colour = 'black',size = 10),
#         axis.text.y = element_text(colour = 'black',size = 10),
#         legend.title = element_text(colour = 'black',size = 10))+
#   labs(x='',y='Fusion sample number',fill='',title = 'Recurrent fusion gene paires')+
#   scale_fill_manual(values = cols)+coord_flip()+ylim(0,20)
# p
# #export::graph2ppt(p,'Fig_S7C_Frequency_of_Fusion_paires_across_clinical_subtype.ppt',height = 5, width = 5,append = F)
# ggsave(p,filename = '/results/Fig_S7C_Frequency_of_Fusion_paires_across_clinical_subtype.pdf',p,width = 5,height = 4)



#### ----------------Fig.S7 D---------------- ####
# df <- apply(Fusion_gene_pair_ProteinCoding,MARGIN = 2,function(x){sum(x)})%>% as.data.frame() 
# colnames(df)[1] <- 'transcripts_number'
# df$Gene_fusion_ID <- rownames(df)
# df <- left_join(meta_fusion_CBCGA,df,by='Gene_fusion_ID')
# df$PAM50_classifier[df$group=='N'] <- 'PT'
# df$Clinical_Subtype[df$group=='N'] <- 'PT'
# 
# df$Clinical_Subtype <- fct_relevel(df$Clinical_Subtype,'HR+HER2-','HR+HER2+','HR-HER2+','TNBC','PT')
# cols <- c('HR+HER2-'='#0085c4','HR+HER2+'='#7ab801','HR-HER2+'='#f2af01','TNBC'='#dc5035','PT' = '#96B9DC' )
# 
# p <- ggboxplot(df, x="Clinical_Subtype", y="transcripts_number", color = "Clinical_Subtype",
#                add = "jitter", add.params = list(size = 0.5, jitter = 0.2))+
#   stat_compare_means(label.y = 65,label.x = 'HR+HER2+')+
#   labs(x='',y='Fusion transcripts per sample',color='')+
#   scale_color_manual(values = cols)+
#   ylim(-1,65)+
#   theme(axis.text.x = element_text(angle = 30, hjust =1, vjust =1,colour = 'black',size = 10))
# p <- ggarrange(p,legend = 'right')
# p
# #export::graph2ppt(p,'Fig_S7D_Fusion_transcripts_across_clinical_subtype.ppt',height = 5, width = 5,append = F)
# ggsave(p,filename = '/results/Fig_S7D_Fusion_transcripts_across_clinical_subtype.pdf',p,width = 5,height = 4)



#### ----------------Fig.S7 E---------------- ####
# cols <- c('LumA'='#0077c9','LumB'='#74d2e8','Her2'='#7552cd','Basal'='#e4002c','Normal'='#cecece','PT' = '#96B9DC')
# df$PAM50_classifier <- fct_relevel(df$PAM50_classifier,c('LumA','LumB','Her2','Basal','Normal','PT'))
# 
# p <- ggboxplot(df, x="PAM50_classifier", y="transcripts_number", color = "PAM50_classifier",
#                add = "jitter", add.params = list(size = 0.5, jitter = 0.2))+
#   stat_compare_means(label.y = 65,label.x = 'LumB')+
#   labs(x='',y='Fusion transcripts per sample',color='')+
#   scale_color_manual(values = cols)+
#   ylim(-1,65)+
#   theme(axis.text.x = element_text(angle = 30, hjust =1, vjust =1,colour = 'black',size = 10))
# p <- ggarrange(p,legend = 'right')
# p
# #export::graph2ppt(p,'Fig_S7E_Fusion_transcripts_across_PAM50.ppt',height = 5, width = 5,append = F)
# ggsave(p,filename = '/results/Fig_S7E_Fusion_transcripts_across_PAM50.pdf',p,width = 5,height = 4)



#### ----------------Fig.S7 F---------------- ####
# chr17q_fusion$reading_frame_arriba[chr17q_fusion$reading_frame_arriba %in% '.'] <- 'Others'
# data <- table(chr17q_fusion$reading_frame_arriba) %>% as.data.frame()
# data$percent <- prop.table(table(chr17q_fusion$reading_frame_arriba))
# data$percent <- round(data$percent,2)
# data$percent <- data$percent*100
# labs <- paste0(data$Var1, " (", data$percent, "%)")
# p <- ggpie(data=data,'Freq',label = labs,lab.pos = "out",
#            fill = 'Var1',color='white',size = 1,
#            palette = c('in-frame'='#1B79B5','Others'='#C9B2D5','out-of-frame'='#E3131E','stop-codon'='#683E95'))
# p
# #export::graph2ppt(p,'Fig_S7F_Fusions_proximal_to_ERBB2.ppt',height = 5, width = 5,append = F)
# ggsave(p,filename = '/results/Fig_S7F_Fusions_proximal_to_ERBB2.pdf',p,width = 5,height = 4)


#### ----------------Fig.S7 G---------------- ####
gf_TT <- Fusion_transcript_ProteinCoding[,TT_id]
gf_TT$exist_in_sample <- apply(gf_TT,1,function(x){sum(x!=0)})
gf_TT$fusion_full <- rownames(gf_TT)

gene <- 'ERBB2'
input_df <- gf_TT[str_detect(gf_TT$fusion_full,gene),c('fusion_full','exist_in_sample')]
df <- str_split(input_df$fusion_full,'[_:]',simplify = T) %>% as.data.frame()
colnames(df) <- c('Gene_pair','Gene1_chr','Gene1_break','Gene1_strand','Gene2_chr','Gene2_break','Gene2_strand','Gene1_ENSG_0','Gene2_ENSG_0')
df$exist_in_sample <- input_df$exist_in_sample
df <- df[!duplicated(df$Gene_pair),]

df1 <- data.frame('chr'=df$Gene1_chr,'start'=as.numeric(df$Gene1_break),'end'=as.numeric(df$Gene1_break),'value'=df$exist_in_sample)
df2 <- data.frame('chr'=df$Gene2_chr,'start'=as.numeric(df$Gene2_break),'end'=as.numeric(df$Gene2_break),'value'=df$exist_in_sample)

circos.clear()
pdf(paste0('/results/FigS7G_',gene,'circos.pdf'),height=5,width=5)
circos.initializeWithIdeogram()
circos.genomicLink(
  df1, df2, border = NA,lwd = 3,
  col = rand_color(nrow(df1),
                   #luminosity = "random", 
                   transparency = 0.5)
)
text(0,0,gene)
graphics.off()

#### ----------------Fig.S7H ---------------- ####
# length(unique(ERBB2_fusion$sample_id))
# length(unique(ERBB2_fusion$fusion_full))
# ERBB2_fusion$kinase <- ifelse(ERBB2_fusion$reading_frame_arriba == 'in-frame' & str_detect(ERBB2_fusion$retained_protein_domains_arriba,'kinase_domain'),'Kinase_in-frame','Others')
# data <- table(ERBB2_fusion$kinase) %>% as.data.frame()
# data$percent <- prop.table(table(ERBB2_fusion$kinase))
# data$percent <- round(data$percent,2)
# data$percent <- data$percent*100
# labs <- paste0(data$Var1, " (", data$percent, "%)")
# p <- ggpie(data=data,'Freq',label = labs,lab.pos = "out",
#            fill = 'Var1',color='white',size = 1,
#            palette = c('Kinase_in-frame'='#1B79B5','Others'='#E3131E'))
# p
# ggsave(p,filename = '/results/Fig_S7H_Fusions_of_ERBB2_pie.pdf',p,width = 5,height = 4)


# #### ----------------Fig.S7I---------------- ####
# ## propensity score matching 
# ERBB2_surv <- read.csv('ERBB2_fusion_PSM_survival.csv',header = T)
# TT_id_ERBB2fusion <- ERBB2_fusion$sample_id
# ERBB2_surv$ERBB2 <- 'No ERBB2 fusion'
# ERBB2_surv$ERBB2[ERBB2_surv$PatientCode %in% meta_fusion_CBCGA$PatientCode[meta_fusion_CBCGA$Gene_fusion_ID %in% TT_id_ERBB2fusion]] <- 'ERBB2 fusion'
# table(ERBB2_surv$ERBB2)
# 
# survdiff(Surv(RFS_months, RFS_status) ~ ERBB2,data=ERBB2_surv)
# p <- ggsurvplot(survfit(Surv(RFS_months, RFS_status) ~ ERBB2,data=ERBB2_surv),
#                 surv.median.line = "hv",
#                 pval = T, 
#                 risk.table = F,
#                 conf.int = FALSE,
#                 legend.title='ERBB2 fusion',
#                 legend.labs=c("ERBB2 fusion","No ERBB2 fusion"),
#                 palette = c("#0077BF","#BE5E55"),
#                 xlim = c(0,120),
#                 break.x.by = 24)+
#   labs(x='Months',y='RFS')
# p
# export::graph2ppt(p$plot,'Fig_S7H_ERBB2_PSM_RFS.ppt',height = 5, width = 5,append = F)
# ggsave(filename = 'Fig_S7H_ERBB2_PSM_RFS.pdf',plot=p$plot,width=5,height=5)



