import pandas as pd
def defineTestSets(df_train,y,whichfeats):
    from sklearn import preprocessing
    df_train_reshuffled = df_train.copy().sample(frac=1, random_state=0).reset_index(drop=True)
    y_reshuffled=y.copy().sample(frac=1, random_state=0).reset_index(drop=True)
    enc = preprocessing.OneHotEncoder()
    X = df_train_reshuffled
    if 'IHC' in whichfeats:
        enc.fit(X.loc[:,['IHC']].values)
        data_x_numeric = enc.transform(X.loc[:,['IHC']].values).toarray()
        data_x_numeric = pd.DataFrame(data_x_numeric,columns=['HR+HER2+','HR+HER2-','TNBC','HR-HER2+'])
        X = X.drop(columns=['IHC'])
        X=pd.concat([X,data_x_numeric],axis=1)
    if 'Clin' in whichfeats:
        enc.fit(X.loc[:,['ajcc']].values)
        data_x_numeric = enc.transform(X.loc[:,['ajcc']].values).toarray()
        data_x_numeric = pd.DataFrame(data_x_numeric,columns=['1','2','3'])
        X = X.drop(columns=['ajcc'])
        X = pd.concat([X,data_x_numeric],axis=1)
    X = X.iloc[:,2:].copy()
    patID = df_train_reshuffled['ID'].copy()
    ytrain=y_reshuffled.iloc[:,2:].copy()
    return X, patID,ytrain

def defineResponse(df):
    df_reshuffled = df.copy().sample(frac=1, random_state=0).reset_index(drop=True)
    ystatus = df_reshuffled['Status'].copy()
    ytime = df_reshuffled['Survival_in_days'].copy()
    return ystatus,ytime

def plotStyle():
    from matplotlib import rcParams
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['.Helvetica Neue DeskInterface']
    rcParams['font.size'] = 18
    rcParams['axes.linewidth'] = 2
    rcParams['grid.linewidth'] = 2
    rcParams['grid.color'] = 'gainsboro'
    rcParams['font.weight'] = 'normal'
    rcParams['axes.labelweight'] = 'bold'
    rcParams['axes.labelsize'] = 21
    rcParams['legend.edgecolor'] = 'none'
    rcParams["axes.spines.right"] = False
    rcParams["axes.spines.top"] = False



def main(whichFeats, rcut,knum,random_seed):
## SETUP
    from CBCGA_validation import  final_test, test_all_models
    import pandas as pd
    import numpy as np
    import pickle
    plotStyle()

    print('Running {} models, with RS={}'.format(whichFeats, random_seed))


    import rpy2.robjects as robjects

    robjects.r.source('D:/zh/ML/trial.R')
    robjects.r.runl(random_seed, whichFeats)
    x = pd.read_csv('D:/zh/ML/inputs/{}Train.csv'.format(whichFeats))
    y = pd.read_csv('D:/zh/ML/inputs/{}yTrain.csv'.format(whichFeats))
## INPUT
    df_test=pd.read_csv('D:/zh/ML/inputs/{}Test.csv'.format(whichFeats))
    df_ytest=pd.read_csv('D:/zh/ML/inputs/{}yTest.csv'.format(whichFeats))



## DATASET
    Xtest_pCR, trialID_pCR,ytest_pCR = defineTestSets(df_test, df_ytest, whichFeats)
    Xfeats=pd.read_csv('D:/zh/ML/CBCGA_mode2/Mode_feats/{}_rcut{}_knum{}_rd{}.csv'.format(whichFeats,rcut,knum,random_seed))

    allfeat=Xfeats.iloc[:,1:].values
    allfeat=allfeat.reshape((allfeat.shape[0]*allfeat.shape[1]))

    ## Limited dataset
    Xnewtest=Xtest_pCR.loc[:,allfeat]


    if 'Clin' in whichFeats:
       Xnewtest = pd.concat([Xnewtest, Xtest_pCR.loc[:, ['1','2','3']]], axis=1)
    if 'IHC' in whichFeats:
       Xnewtest = pd.concat([Xnewtest, Xtest_pCR.loc[:, ['HR+HER2+','HR+HER2-','TNBC','HR-HER2+']]],axis=1)
## MODELS
    ### pCR
    stamp = 'submission'
    outdir = 'results_{}_r{}_k{}_rs{}'.format(stamp, rcut, knum,random_seed)
    pcr_refits = pickle.load(open('D:/zh/ML/CBCGA_mode2/datav/{}/{}_pcr_models.pkl'.format(outdir, whichFeats), 'rb'))

    y_catarray = ytest_pCR['Status']
    y_cattime = ytest_pCR['Survival_in_days']
    ynew = np.dtype({'names': ['Status', 'Survival_in_days'], 'formats': ['?', '<f2']})
    yrun = np.array(list(zip(y_catarray, y_cattime)), dtype=ynew)

    test_all_models(Xnewtest, yrun, pcr_refits, 'pCR_{}_r{}'.format(whichFeats, rcut),rcut,random_seed,whichFeats,knum)

if __name__ == "__main__":

    main('RNA_Path_IHC_Clin',0.6,5,111)