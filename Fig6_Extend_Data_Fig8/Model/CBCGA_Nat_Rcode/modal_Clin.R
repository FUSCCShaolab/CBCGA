#Singlemodal_Clin
setwd('C:/Users/lenovo/PycharmProjects/pythonProject')
CBCGA_Cohort.Info=read.csv('extended.csv',header = T)
intersects <- function (...) {
  Reduce(intersect, list(...))
}
#需要变动不同模态
Trainset_IHC=outputlist$Train_IHC
Trainset_Clin=outputlist$Train_Clinc
Trainset_Met=outputlist$Train_met
Trainset_RNA=outputlist$Train_RNA
Trainset_rad=outputlist$Train_radio
Trainset_Path=outputlist$Train_path
test_IHC=outputlist$Test_IHC
test_Clin=outputlist$Test_Clinic
test_Met=outputlist$Test_met
test_rad=outputlist$Test_radio
test_path=outputlist$Test_path
test_rna=outputlist$Test_RNA
#只需改动
id=intersects(rownames(Trainset_Clin))
Train=as.data.frame(cbind(Trainset_Clin[match(id,rownames(Trainset_Clin)),]))
rownames(Train)=id
colnames(Train)='ajcc'
Test=as.data.frame(cbind(test_Clin))
rownames(Test)=rownames(test_Clin)
colnames(Test)='ajcc'

#修改不合适列名 不需要修改
if(is.element('Stage (pN)',colnames(Train))){
  colnames(Train)[match('Stage (pN)',colnames(Train))]='pN'
}
if(is.element('Stage (pT)',colnames(Train))){
  colnames(Train)[match('Stage (pT)',colnames(Train))]='pT'
}
if(is.element('Trainset_IHC[match(id, rownames(Trainset_IHC)), ]',colnames(Train))){
  colnames(Train)[match('Trainset_IHC[match(id, rownames(Trainset_IHC)), ]',colnames(Train))]='IHC'
}
if(is.element('Stage (pN)',colnames(Test))){
  colnames(Test)[match('Stage (pN)',colnames(Test))]='pN'
}
if(is.element('Stage (pT)',colnames(Test))){
  colnames(Test)[match('Stage (pT)',colnames(Test))]='pT'
}
if(is.element('Trainset_IHC[match(id, rownames(Trainset_IHC)), ]',colnames(Test))){
  colnames(Test)[match('Trainset_IHC[match(id, rownames(Trainset_IHC)), ]',colnames(Test))]='IHC'
}
if(is.element("Trainset_Clin[match(id, rownames(Trainset_Clin)), ]",colnames(Train))){
  colnames(Train)[match("Trainset_Clin[match(id, rownames(Trainset_Clin)), ]",colnames(Train))]='ajcc'
}
#不需要修改
Trainrfs=CBCGA_Cohort.Info$RFS_status[match(id,CBCGA_Cohort.Info$PatientCode)]
Trainrfstime=CBCGA_Cohort.Info$RFS_time_Months[match(id,CBCGA_Cohort.Info$PatientCode)]
Testrfs=CBCGA_Cohort.Info$RFS_status[match(rownames(Test),CBCGA_Cohort.Info$PatientCode)]
Testrfstime=CBCGA_Cohort.Info$RFS_time_Months[match(rownames(Test),CBCGA_Cohort.Info$PatientCode)]

yTrain=cbind(Trainrfstime,Trainrfs)
colnames(yTrain)=c('Survival_in_days','Status')
yTest=cbind(Testrfstime,Testrfs)
colnames(yTest)=c('Survival_in_days','Status')