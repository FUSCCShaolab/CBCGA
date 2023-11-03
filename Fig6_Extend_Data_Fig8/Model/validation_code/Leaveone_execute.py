
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
        enc.fit(X.loc[:, ['IHC']].values)
        data_x_numeric = enc.transform(X.loc[:, ['IHC']].values).toarray()
        data_x_numeric = pd.DataFrame(data_x_numeric, columns=['HR+HER2+', 'HR+HER2-', 'TNBC', 'HR-HER2+'])
        X = X.drop(columns=['IHC'])
        X = pd.concat([X, data_x_numeric], axis=1)
    X = X.iloc[:,:].copy()
    patID = y_reshuffled['ID'].copy()
    ytrain=y_reshuffled.iloc[:,2:].copy()
    return X, patID,ytrain

def main(whichFeats, rcut,knum,random_seed):
    from CBCGA_validation import final_test, test_all_models
    import pandas as pd
    import pickle
    print('Running {} models, with RS={}'.format(whichFeats, random_seed))
    x = pd.read_csv('/Users/zhanghang/Desktop/test80_geneset/RNA_Met_Path_IHC_ClinTest.csv')
    df_ytest = pd.read_csv('/Users/zhanghang/Desktop/test80_geneset/RNA_Met_Path_IHC_ClinyTest.csv')
    Feats = ['M168T190.POS', 'M297T233.POS', 'M149T91.NEG', 'M773T388.NEG', 'M778T361.NEG',
             'T.GLCM.IntensityMax.Mean', 'T.GLCM.IntensityMax.Skew',
             'T.Morph.CellEccentricities.Kurtosis', 'T.GLCM.IntensityMean.Kurtosis', 'Normal', 'GFRA1', 'TFF3',
             'CCNB1', 'CTSV', 'CDH3',
             'Clin', 'IHC']
    loo_feats = [x for x in Feats if x != whichFeats]
    df_test = x.loc[:, loo_feats]
    df_testnew, df_id, df_ytrain = defineTrainingSets(df_test, df_ytest, loo_feats)
    import numpy as np
    y_catarray = df_ytrain['Status']
    y_cattime = df_ytrain['Survival_in_days']
    ynew = np.dtype({'names': ['Status', 'Survival_in_days'], 'formats': ['?', '<f2']})
    yrun = np.array(list(zip(y_catarray, y_cattime)), dtype=ynew)
    stamp = 'submission'
    outdir = 'results_{}_r{}_k{}_rs{}'.format(stamp, rcut, knum, random_seed)
    pcr_refits = pickle.load(open('/Users/zhanghang/PycharmProjects/pythonProject/CBCGA_mode2/leaveone/{}/{}_pcr_refits.pkl'.format(outdir,whichFeats),'rb'))
    print(df_testnew)
    test_all_models(df_testnew, yrun, pcr_refits, 'pCR_{}_r{}'.format(whichFeats, rcut), rcut, random_seed, whichFeats,knum, df_id)


if __name__ == "__main__":
    main('ZH', 0.6, 5, 111)