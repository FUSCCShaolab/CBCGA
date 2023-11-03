rm(list = ls()) ; graphics.off()
#------------------------------------------------------------------------------------#
# package prepare
#------------------------------------------------------------------------------------#
library(coin)
library(ggrepel)
library(plyr)
library(tidyverse)
library(export)
library(pheatmap)
load("/data/CBCGA.Extended_MergedData_V2.5_220722.Rdata")
load("/data/Fig3_S3_data.Rdata")

#--------------------------------------------------------------------
# part1: the whole population
#--------------------------------------------------------------------
########### data input -------
##### CBCGA
thre = CBCGA_GISTICgene.thre
gene.loc = gene.loc.noX
thre = thre[rownames(gene.loc),]

clinic_CBCGA = CBCGA_Cohort.Info
clinic_CBCGA = clinic_CBCGA[intersect(rownames(clinic_CBCGA), colnames(thre)),]
clinic_CBCGA = clinic_CBCGA[clinic_CBCGA$`Histological type` == "IDC",]

subtype_CBCGA = list(HRposHER2neg = rownames(clinic_CBCGA)[(clinic_CBCGA$`ER status` == "Positive" | clinic_CBCGA$`PR status` == "Positive") & clinic_CBCGA$`HER2 status (combined)` == "Negative"],
                     HRposHER2pos = rownames(clinic_CBCGA)[(clinic_CBCGA$`ER status` == "Positive" | clinic_CBCGA$`PR status` == "Positive") & clinic_CBCGA$`HER2 status (combined)` == "Positive"],
                     HRnegHER2pos = rownames(clinic_CBCGA)[(clinic_CBCGA$`ER status` == "Negative" & clinic_CBCGA$`PR status` == "Negative") & clinic_CBCGA$`HER2 status (combined)` == "Positive"],
                     HRnegHER2neg = rownames(clinic_CBCGA)[ clinic_CBCGA$`ER status` == "Negative" & clinic_CBCGA$`PR status` == "Negative" & clinic_CBCGA$`HER2 status (combined)` == "Negative"]
)

##### TCGA
clinic_TCGA = clinic_TCGA[clinic_TCGA$Race.Category == "WHITE" & !is.na(clinic_TCGA$Race.Category), ]
clinic_TCGA = clinic_TCGA[substring(rownames(clinic_TCGA),14,15) == "01",]
rownames(clinic_TCGA) = substring(rownames(clinic_TCGA),1,12)
clinic_TCGA = clinic_TCGA[clinic_TCGA$histological_type == "Infiltrating Ductal Carcinoma", ]

tmp = intersect(rownames(TCGA_cnv),rownames(thre))
TCGA_cnv = TCGA_cnv[tmp,]
thre = thre[tmp,]

tmp = intersect(rownames(clinic_TCGA),colnames(TCGA_cnv))
TCGA_cnv = TCGA_cnv[,tmp]
clinic_TCGA = clinic_TCGA[tmp,]

subtype_TCGA = list(HRposHER2neg = na.omit(rownames(clinic_TCGA)[(clinic_TCGA$ER_Final == "Positive" | clinic_TCGA$PR_Final == "Positive") & clinic_TCGA$HER2_Final == "Negative"]),
                    HRposHER2pos = na.omit(rownames(clinic_TCGA)[(clinic_TCGA$ER_Final == "Positive" | clinic_TCGA$PR_Final == "Positive") & clinic_TCGA$HER2_Final == "Positive"]),
                    HRnegHER2pos = na.omit(rownames(clinic_TCGA)[(clinic_TCGA$ER_Final == "Negative" & clinic_TCGA$PR_Final == "Negative") & clinic_TCGA$HER2_Final == "Positive"]),
                    HRnegHER2neg = na.omit(rownames(clinic_TCGA)[(clinic_TCGA$ER_Final == "Negative" & clinic_TCGA$PR_Final == "Negative") & clinic_TCGA$HER2_Final == "Negative"])
)


########### DCG & draw ---------
CBCGA = rownames(clinic_CBCGA)
TCGA = rownames(clinic_TCGA)

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

##### compare ------
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

## gain DCG -------
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
#write.csv(CBCGAVStcga.gain.res,file = paste0("/results/CBCGAVStcga_WHITE_gain_res.csv"))

## loss DCG -------
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
#write.csv(CBCGAVStcga.loss.res,file = paste0("/results/CBCGAVStcga_WHITE_loss_res.csv"))

## gain plot -------
#CBCGAVStcga.gain.res = read.csv(paste0("/results/CBCGAVStcga_WHITE_gain_res.csv"),row.names = 1)
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
                                               (CBCGAVStcga.gain.res$CBCGA > 0.1 & CBCGAVStcga.gain.res$TCGA >0.1),] ,
                  aes(label=rownames(CBCGAVStcga.gain.res[CBCGAVStcga.gain.res$adj.p < 0.05 &
                                                            (CBCGAVStcga.gain.res$CBCGA > 0.1 & CBCGAVStcga.gain.res$TCGA >0.1),])),
                  col="black",alpha = 1,size=3) 
ggsave(p,filename = paste0("/results/CBCGAvsTCGA_WHITE_gain.pdf"),height = 4,width = 4)
#export::graph2ppt(p,file = paste0("/results/CBCGAvsTCGA_WHITE_gain.pdf"),height = 4,width = 4)

p = ggplot(CBCGAVStcga.gain.res,aes(x=CBCGA, y = TCGA)) + 
  geom_point(alpha = 1, size = 2.5, col="#cecece") + 
  geom_point(data = CBCGAVStcga.gain.res[CBCGAVStcga.gain.res$adj.p < 0.05,] , 
             alpha = 1, size = 2.5, col="#C0052A") +
  geom_abline(intercept=0,slope=1, linetype="dashed", size = 1.2, color = '#747475')+
  xlim(0,0.25) + ylim(0,0.25) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1.2)) +
  geom_text_repel(data= CBCGAVStcga.gain.res[CBCGAVStcga.gain.res$adj.p < 0.05 & 
                                               (CBCGAVStcga.gain.res$CBCGA > 0.1 | CBCGAVStcga.gain.res$TCGA > 0.1),] ,
                  aes(label=rownames(CBCGAVStcga.gain.res[CBCGAVStcga.gain.res$adj.p < 0.05 &
                                                            (CBCGAVStcga.gain.res$CBCGA > 0.1 | CBCGAVStcga.gain.res$TCGA > 0.1),])),
                  col="black",alpha = 1,size=3) 
ggsave(p,filename = paste0("/results/CBCGAvsTCGA_WHITE_gain2.pdf"),height = 4,width = 4)
#export::graph2ppt(p,file = paste0("/results/CBCGAvsTCGA_WHITE_gain2.pdf"),height = 4,width = 4)

## loss plot-------
#CBCGAVStcga.loss.res = read.csv(paste0("/results/CBCGAVStcga_WHITE_loss_res.csv"),row.names = 1)
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
                                               (CBCGAVStcga.loss.res$CBCGA > 0.1 | CBCGAVStcga.loss.res$TCGA > 0.1 ),] ,
                  aes(label=rownames(CBCGAVStcga.loss.res[CBCGAVStcga.loss.res$adj.p < 0.05 &
                                                            (CBCGAVStcga.loss.res$CBCGA > 0.1 | CBCGAVStcga.loss.res$TCGA > 0.1 ),])),
                  col="black",alpha = 1,size=3,
                  max.overlaps = 20) 
ggsave(p,filename = paste0("/results/CBCGAvsTCGA_WHITE_loss.pdf"),height = 4,width = 4)
#export::graph2ppt(p,file = paste0("/results/CBCGAvsTCGA_WHITE_loss.pdf") ,height = 4,width = 4)

p = ggplot(CBCGAVStcga.loss.res,aes(x=CBCGA, y = TCGA)) + 
  geom_point(alpha = 1, size = 2.5, col="#cecece") + 
  geom_point(data = CBCGAVStcga.loss.res[CBCGAVStcga.loss.res$adj.p < 0.05,] , 
             alpha = 1, size = 2.5, col="#4E85AC") +
  geom_abline(intercept=0,slope=1, linetype="dashed", size = 1.2, color = '#747475')+
  xlim(0,0.25) + ylim(0,0.25) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1.2)) +
  geom_text_repel(data= CBCGAVStcga.loss.res[CBCGAVStcga.loss.res$adj.p < 0.05 & 
                                               (CBCGAVStcga.loss.res$CBCGA > 0.1 | CBCGAVStcga.loss.res$TCGA > 0.1 ),] ,
                  aes(label=rownames(CBCGAVStcga.loss.res[CBCGAVStcga.loss.res$adj.p < 0.05 &
                                                            (CBCGAVStcga.loss.res$CBCGA > 0.1 | CBCGAVStcga.loss.res$TCGA > 0.1 ),])),
                  col="black",alpha = 1,size=3,
                  max.overlaps = 20) 
ggsave(p,filename = paste0("/results/CBCGAvsTCGA_WHITE_loss2.pdf"),height = 4,width = 4)
#export::graph2ppt(p,file = paste0("/results/CBCGAvsTCGA_WHITE_loss2.pdf") ,height = 4,width = 4)

#--------------------------------------------------------------------
# part1: IHC subtype
#--------------------------------------------------------------------
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
  ggsave(p,filename = paste0("/results/CBCGAvsTCGA_WHITE_gain_",k,".pdf"),height = 4,width = 4)
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
  ggsave(p,filename = paste0("/results/CBCGAvsTCGA_WHITE_loss_",k,".pdf"),height = 4,width = 4)
  #export::graph2ppt(p,file = paste0("/results/CBCGAvsTCGA_WHITE_loss_",k,".pdf"),height = 4,width = 4)
}


