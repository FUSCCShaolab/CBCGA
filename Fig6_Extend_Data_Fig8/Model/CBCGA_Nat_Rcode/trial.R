#CBCGA_Natmodel_run
runl=function(rd,mod){
  setwd('D:/zh/ML')
  CBCGA_Cohort.Info=read.csv('CBCGA_Cohort_info.csv',header = T,row.names = 1)
  z1=read.csv('z1.csv',header = T,row.names = 1)
  set.seed(rd)
  sub_test=sampling::strata(z1,stratanames ='RFS', size=c(66,14), method = 'srswor')
  data_test=z1[sub_test$ID_unit,]
  data_train=z1[-sub_test$ID_unit,]
  out=list()
  setwd('D:/zh/ML')
  pathologynew=read.csv('CBCGA_pathologynew.csv',header = T,row.names = 1)
  radiomicnew=read.csv('CBCGA_radiomicnew.csv',header = T,row.names = 1)
  rnanew=read.csv('CBCGA_rnanew.csv',header = T,row.names = 1)
  metabolismnew=read.csv('CBCGA_metabolismnew.csv',header = T,row.names = 1)
  

  a1=c()
  for (k in 1:length(data_test$ID)) {
    a1=append(a1,which(rownames(pathologynew)==data_test$ID[k])) 
  }
  Trainset_candidate_Path=pathologynew[-a1,]
  Test_can_patho=pathologynew[a1,]
  a2=c()
  for (k in 1:length(data_test$ID)) {
    a2=append(a2,which(rownames(radiomicnew)==data_test$ID[k])) 
  }
  Trainset_candidate_rad=radiomicnew[-a2,]
  Test_can_rad=radiomicnew[a2,]
  a3=c()
  for (k in 1:length(data_test$ID)) {
    a3=append(a3,which(rownames(rnanew)==data_test$ID[k])) 
  }
  Trainset_candidate_RNA=rnanew[-a3,]
  Test_can_RNA=rnanew[a3,]
  a4=c()
  for (k in 1:length(data_test$ID)) {
    a4=append(a4,which(rownames(metabolismnew)==data_test$ID[k])) 
  }
  Trainset_candidate_Met=metabolismnew[-a4,]
  Test_can_Met=metabolismnew[a4,]
  a5=c()
  for (k in 1:length(data_test$ID)) {
    a5=append(a5,which(CBCGA_Cohort.Info$PatientCode==data_test$ID[k])) 
  }
  Trainset_candidate_Clin=as.data.frame(CBCGA_Cohort.Info[-a5,36])
  Test_can_Clin=as.data.frame(CBCGA_Cohort.Info[a5,36])
  rownames(Trainset_candidate_Clin)=CBCGA_Cohort.Info$PatientCode[-a5]
  colnames(Trainset_candidate_Clin)='ajcc'
  rownames(Test_can_Clin)=CBCGA_Cohort.Info$PatientCode[a5]
  colnames(Test_can_Clin)='ajcc'
  a6=c()
  for (k in 1:length(data_test$ID)) {
    a6=append(a6,which(CBCGA_Cohort.Info$PatientCode==data_test$ID[k])) 
  }
  
  Trainset_candidate_IHC=as.data.frame(CBCGA_Cohort.Info[-a6,46])
  Test_can_IHC=as.data.frame(CBCGA_Cohort.Info[a6,46])
  rownames(Trainset_candidate_IHC)=rownames(CBCGA_Cohort.Info[-a6,])
  colnames(Trainset_candidate_IHC)='IHC'
  colnames(Test_can_IHC)='IHC'
  outputlist<<-list(Train_met=Trainset_candidate_Met,Train_path=Trainset_candidate_Path,Train_RNA=Trainset_candidate_RNA,
                    Train_radio=Trainset_candidate_rad,Train_Clinc=Trainset_candidate_Clin,Train_IHC=Trainset_candidate_IHC,
                    Test_met=Test_can_Met,Test_path=Test_can_patho,Test_RNA=Test_can_RNA,Test_radio=Test_can_rad,
                    Test_Clinic=Test_can_Clin,Test_IHC=Test_can_IHC)
  setwd('D:/zh/ML')
  source(paste('D:/zh/ML/CBCGA_Nat_Rcode/modal_',mod,'.R',sep = ''))
  idtrain=rownames(Train)
  Train=cbind(ID=idtrain,Train)
  yTrain=as.data.frame(yTrain)
  yTrain=cbind(ID=idtrain,yTrain)
  idtest=rownames(data_test)
  Test=cbind(ID=idtest,Test)
  yTest=as.data.frame(yTest)
  yTest=cbind(ID=idtest,yTest)
  rownames(Train)=c(1:length(idtrain))
  rownames(yTrain)=c(1:length(idtrain))
  rownames(Test)=c(1:dim(Test)[1])
  rownames(yTest)=c(1:dim(Test)[1])
  setwd('D:/zh/ML/inputs')
  write.csv(Train,paste(mod,'Train.csv',sep = ''))
  write.csv(yTrain,paste(mod,'yTrain.csv',sep = ''))
  write.csv(Test,paste(mod,'Test.csv',sep = ''))
  write.csv(yTest,paste(mod,'yTest.csv',sep = ''))
  return()
  
}
