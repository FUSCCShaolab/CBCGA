import pandas as pd
def defineSplits(X,ycateg,random_state):
    from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=int(random_state))
    splits = []
    for (tr,ts) in cv.split(X, ycateg):
        splits.append((tr,ts))
    return splits

def defineTrainingSets(df_train,y,whichfeats):
    import pandas as pd
    from sklearn import preprocessing
    df_train_reshuffled = df_train.copy().sample(frac=1, random_state=0).reset_index(drop=True)
    y_reshuffled=y.copy().sample(frac=1, random_state=0).reset_index(drop=True)
    enc = preprocessing.OneHotEncoder()
    X = df_train_reshuffled
    if 'Clin' in whichfeats:
        enc.fit(X.loc[:,['Clin']].values)
        data_x_numeric = enc.transform(X.loc[:,['Clin']].values).toarray()
        data_x_numeric = pd.DataFrame(data_x_numeric,columns=['1','2','3'])
        X = X.drop(columns=['Clin'])
        X = pd.concat([X,data_x_numeric],axis=1)
    if 'IHC' in whichfeats:
        enc.fit(X.loc[:,['IHC']].values)
        data_x_numeric = enc.transform(X.loc[:,['IHC']].values).toarray()
        data_x_numeric = pd.DataFrame(data_x_numeric,columns=['HR+HER2+','HR+HER2-','TNBC','HR-HER2+'])
        X = X.drop(columns=['IHC'])
        X=pd.concat([X,data_x_numeric],axis=1)
    X = X.iloc[:,:].copy()
    patID = y_reshuffled['ID'].copy()
    ytrain=y_reshuffled.iloc[:,2:].copy()
    return X, patID,ytrain

def main(whichFeats, rcut, knum,random_state):
       from CBCGA_classification import run_all_models,refit_all_models
       import pickle
       allfeatcsv = pd.read_csv('/Users/zhanghang/Desktop/test80_geneset/RNA_Met_Path_IHC_ClinTrain.csv')
       y = pd.read_csv('/Users/zhanghang/Desktop/test80_geneset/RNA_Met_Path_IHC_ClinyTrain.csv')
       Feats = ['M168T190.POS', 'M297T233.POS', 'M149T91.NEG', 'M773T388.NEG', 'M778T361.NEG',
                'T.GLCM.IntensityMax.Mean', 'T.GLCM.IntensityMax.Skew',
                'T.Morph.CellEccentricities.Kurtosis', 'T.GLCM.IntensityMean.Kurtosis', 'Normal', 'GFRA1', 'TFF3',
                'CCNB1', 'CTSV', 'CDH3',
                'Clin', 'IHC']
       loo_feats = [x for x in Feats if x != whichFeats]
       x = allfeatcsv.loc[:, loo_feats]
       df_train, df_id, df_ytrain = defineTrainingSets(x, y, loo_feats)
       print(df_train)
       import numpy as np
       y_catarray = df_ytrain['Status']
       y_cattime = df_ytrain['Survival_in_days']
       ynew = np.dtype({'names': ['Status', 'Survival_in_days'], 'formats': ['?', '<f2']})
       yrun = np.array(list(zip(y_catarray, y_cattime)), dtype=ynew)
       splits = defineSplits(df_train, df_ytrain['Status'], random_state)
       print('spllit done')
       pcr_models = run_all_models(df_train, yrun, splits, float(rcut),rd=random_state)  # 分别是X y splits=cv cut=cut cut是dropcollinear中的corr cut
       pcr_refits = refit_all_models(df_train, yrun, pcr_models, splits, whichFeats,'c_index', rcut,random_state, df_id)
       import os
       stamp = 'submission'
       outdir = 'results_{}_r{}_k{}_rs{}'.format(stamp, rcut, knum, random_state)
       if not os.path.exists('/Users/zhanghang/PycharmProjects/pythonProject/CBCGA_mode2/leaveone/{}'.format(outdir)):
              os.makedirs('/Users/zhanghang/PycharmProjects/pythonProject/CBCGA_mode2/leaveone/{}'.format(outdir))
              print('Making: /Users/zhanghang/PycharmProjects/pythonProject/CBCGA_mode2/leaveone/{}'.format(outdir))
       with open('/Users/zhanghang/PycharmProjects/pythonProject/CBCGA_mode2/leaveone/{}/{}_pcr_models.pkl'.format(outdir,whichFeats),'wb') as f:
              pickle.dump(pcr_models, f)
       with open('/Users/zhanghang/PycharmProjects/pythonProject/CBCGA_mode2/leaveone/{}/{}_pcr_refits.pkl'.format(outdir,whichFeats),'wb') as f:
              pickle.dump(pcr_refits, f)
       f.close()

if __name__ == "__main__":
    main('ZH',0.6,5,111)
