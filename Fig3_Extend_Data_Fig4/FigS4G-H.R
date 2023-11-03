# load SNF cluster annotation in Fig3E

# multivariate Cox regression

library(tableone) 
library(forestplot)
library(stringr)

data <- as.data.frame(SNF_Cluster)

# rename according to figure
data[data$SNF_Cluster==1,"Multi_omics_cluster"] <- "cluster2" 
data[data$SNF_Cluster==2,"Multi_omics_cluster"] <- "cluster4" 
data[data$SNF_Cluster==4,"Multi_omics_cluster"] <- "cluster1" 
data[data$SNF_Cluster==3,"Multi_omics_cluster"] <- "cluster3" 

data$PatientCode <- rownames(data)
data <- merge(data,CBCGA_Cohort.Info[,c(1,10,31,33,44,52,53)],all.x = T,by = "PatientCode")
colnames(data)[4:9] <- c("PAM50","size","LN","Clinical_Subtype","RFS_status","RFS_months")

# data <- data[data$Clinical_Subtype=="HR+HER2-",]

data[which(data$size<=2),"size"] <- 0
data[which(data$size>2),"size"] <- 1
data$size <- as.factor(data$size)

data[which(data$LN==0),"LN"] <- 0
data[which(data$LN!=0),"LN"] <- 1
data$LN <- as.factor(data$LN)

data$PAM50 <- as.factor(data$PAM50)
data$PAM50 <- relevel(data$PAM50,ref = "LumA")

for (i in 8:9) {
  data[,i] <- as.numeric(data[,i])
}

res.cox_RFS <- coxph(Surv(RFS_months, RFS_status) ~ PAM50 + size + LN, data = data)

# forest plot
mul_cox1 <- summary(res.cox_RFS)
multi1<-as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))
multi2<-ShowRegTable(res.cox_RFS,exp=TRUE,digits=2,pDigits =3,printToggle = TRUE, quote=FALSE, ciFun=confint)

result <-cbind(multi1,multi2)
result<-tibble::rownames_to_column(result, var = "Characteristics")

ins <- function(x) {c(x, rep(NA,ncol(result)-1))}

for(i in 5:6) {
  result[,i] <- as.character(result[,i])
}

result<-rbind(c("Characteristics", NA, NA, NA, "HR(95%CI)","p"),
              ins("PAM50"),
              ins("LumA"),  
              result[c(3,2,1,4),], 
              ins("Size"),
              ins("¡Ü2cm"),
              result[5,],
              ins("Lymph Node Status"),
              ins("Negative"),
              result[6,],
              c(NA, NA, NA, NA, NA,NA))

for(i in 2:4) {
  result[,i] <- as.numeric(result[,i])
}

myVars <- c("PAM50","Size","Lymph Node Status")
catVars <-  c("PAM50","Size","Lymph Node Status")

colnames(data)[4:6] <- c("PAM50","Size","Lymph Node Status")
table1<- print(CreateTableOne(vars = myVars,
                              data = data[,c(4:6)],
                              factorVars = catVars),
               showAllLevels=TRUE)

N<-rbind(c(NA,NA),table1[c(1,2,5,4,3,6),],c(NA,NA),table1[7:8,],c(NA,NA),table1[9:10,],c(NA,NA))       
N<-as.data.frame(N[,-1])

result1<-cbind(result,N)
result1<-result1[,c(1,7,2:6)]

for(i in 2:7) {
  result1[,i] <-  as.character(result1[,i])
}

result1<-rbind(c("Characteristics","Number(%)",NA,NA,NA,"HR (95%CI)","P.value"),result1[2:nrow(result1),])

for(i in 3:5) {
  result1[,i] <- as.numeric(result1[,i])
}

pdf(file="Fig4G.pdf",width = 15,height = 6)
# pdf(file="Fig4H.pdf",width = 15,height = 6)
forestplot(result1[,c(1,2,6,7)],mean=result1[,3],lower=result1[,4],upper=result1[,5],zero=1,boxsize=0.2,graph.pos= "right",
           hrzl_lines=list("1" = gpar(lty=1,lwd=2),"2" = gpar(lty=2),"13"=  gpar(lwd=2,lty=1,columns=c(1:4))),
           graphwidth = unit(.25,"npc"), xlab="better survival  worse survival",xticks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13), 
           is.summary=c(T,T,F,F,F,F,F,T,F,F,T,F,F), 
           txt_gp=fpTxtGp(label=gpar(cex=1),ticks=gpar(cex=1), xlab=gpar(cex=1),title=gpar(cex=2)),lwd.zero=2,
           lwd.ci=2,lwd.xaxis=2,lty.ci=1,ci.vertices =T,ci.vertices.height=0.1,clip=c(0,13),lineheight=unit(16,'mm'),   
           line.margin=unit(16,'mm'),colgap=unit(6,'mm'),fn.ci_norm="fpDrawDiamondCI",title=NULL,
           col=fpColors(box ='#021eaa',lines ='#021eaa',zero = "black"))
dev.off()
