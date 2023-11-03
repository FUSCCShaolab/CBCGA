# 直接用R来获取输入数据，然后使用 python 交互
# 输入whichfeats:就是模型的名字，然后使用与R交互确定训练集和验证集 注意导出
# X(matrix)
# y(time,status)
# dropcoll--rcut
# randome state
# main(model,rcut,rd)
import pandas as pd
import numpy as np
from sklearn.feature_selection import SelectKBest
from sklearn.preprocessing import StandardScaler
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sklearn.base import BaseEstimator, TransformerMixin, RegressorMixin, clone
from scipy.stats import spearmanr

def defineFeatures(whichFeats):


    metab_feature = pd.read_csv('D:/zh/ML/CBCGA_metabolismnew.csv').columns[1:]
    patho_feature = pd.read_csv('D:/zh/ML/CBCGA_pathologynew.csv').columns[1:]
    radio_feature = pd.read_csv('D:/zh/ML/CBCGA_radiomicnew.csv').columns[1:]
    rna_feature = pd.read_csv('D:/zh/ML/CBCGA_rnanew.csv').columns[1:]

    feature = {}
    if 'Met' in whichFeats:
        feature['Met']=metab_feature
    if 'Path' in whichFeats:
        feature['Path']=patho_feature
    if 'Rad' in whichFeats:
        feature['Rad']=radio_feature
    if 'RNA' in whichFeats:
        feature['RNA']=rna_feature
    return feature


def defineSplits(X,ycateg,random_state):
    from sklearn.model_selection import StratifiedKFold
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=int(random_state))
    splits = []
    for (tr,ts) in cv.split(X, ycateg):
        splits.append((tr,ts))
    return splits

def defineTrainingSets(df_train,y,whichfeats):
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

def fit_and_score_features(X, y):
    n_features = X.shape[1]
    scores = np.empty(n_features)
    m = CoxPHSurvivalAnalysis()
    for j in range(n_features):
        Xj = X[:, j:j+1]
        m.fit(Xj, y)
        scores[j] = m.score(Xj, y)
    return scores
class DropCollinear(BaseEstimator, TransformerMixin):
    def __init__(self, thresh):
        self.uncorr_columns = None
        self.thresh = thresh

    def fit(self, X, y):
        import numpy
        cols_to_drop = []

        # Find variables to remove
        X_corr = X.corr()
        large_corrs = X_corr>self.thresh
        indices = np.argwhere(large_corrs.values)
        indices_nodiag = np.array([[m,n] for [m,n] in indices if m!=n])
        if indices_nodiag.size>0:
            coxm = CoxPHSurvivalAnalysis()
            coxn = CoxPHSurvivalAnalysis()
            indices_nodiag_lowfirst = np.sort(indices_nodiag, axis=1)
            correlated_pairs = np.unique(indices_nodiag_lowfirst, axis=0)
            resp_corrs=np.empty(shape=(0,2))
            for [m,n] in correlated_pairs:
                coxm.fit(X.iloc[:,m].values.reshape(-1,1), y)
                x1=coxm.score(X.iloc[:,m].values.reshape(-1,1), y)
                coxn.fit(X.iloc[:,n].values.reshape(-1,1), y)
                x2=coxn.score(X.iloc[:,n].values.reshape(-1,1), y)
                resp_corrs=numpy.concatenate((resp_corrs,[[x1,x2]]),axis=0)
            element_to_drop = np.argmin(resp_corrs, axis=1)
            list_to_drop = np.unique(correlated_pairs[range(element_to_drop.shape[0]),element_to_drop])
            cols_to_drop = X.columns.values[list_to_drop]

        cols_to_keep = [c for c in X.columns.values if c not in cols_to_drop]
        self.uncorr_columns = cols_to_keep

        return self

    def transform(self, X):
        return X[self.uncorr_columns]

    def get_params(self, deep=False):
        return {'thresh': self.thresh}


class SelectAtMostKBest(SelectKBest):
    def _check_params(self, X, y):
        if not (self.k == "all" or 0 <= self.k <= X.shape[1]):
            self.k = "all"
def main(whichFeats, rcut, knum,random_state):
## SETUP
    from CBCGA_classification import run_all_models, refit_all_models
    import pandas as pd
    import pickle
    import rpy2.robjects as robjects
    print('Running {} models'.format(whichFeats))
    robjects.r.source('D:/zh/ML/trial.R')
    robjects.r.runl(random_state,whichFeats)
    x=pd.read_csv('D:/zh/ML/inputs/{}Train.csv'.format(whichFeats))
    y=pd.read_csv('D:/zh/ML/inputs/{}yTrain.csv'.format(whichFeats))
    df_train, df_id, df_ytrain=defineTrainingSets(x,y,whichFeats)
    import numpy as np
    y_catarray = df_ytrain['Status']
    y_cattime = df_ytrain['Survival_in_days']
    ynew = np.dtype({'names': ['Status', 'Survival_in_days'], 'formats': ['?', '<f2']})
    yrun = np.array(list(zip(y_catarray, y_cattime)), dtype=ynew)
    features=defineFeatures(whichFeats)
    dropcoll=DropCollinear(rcut)
    select_feats={}
    all_col_names=[]
    for feat in features.keys():
        DropedDat=dropcoll.fit_transform(df_train.loc[:,features[feat]],yrun)
        KB=SelectAtMostKBest(fit_and_score_features,k=knum)
        threshold = sorted(KB.scores_,reverse=True)[knum-1]
        col_names=[]
        for i in zip(KB.scores_,KB.feature_names_in_):
            if i[0] >= threshold:
                col_names.append(i[1])
                all_col_names.append(i[1])

    select_feats[whichFeats]=all_col_names
    Xfeats=pd.DataFrame(select_feats)
    Xfeats.to_csv('D:/zh/ML/CBCGA_mode2/Mode_feats/{}_rcut{}_knum{}_rd{}.csv'.format(whichFeats,rcut,knum,random_state))
    allfeat=Xfeats.values
    allfeat=allfeat.reshape((allfeat.shape[0]*allfeat.shape[1]))
    df_selecttrain=df_train.loc[:,allfeat]
    print('r done')
## INPUT

## DATASET
    if 'Clin' in whichFeats:
        df_selecttrain = pd.concat([df_selecttrain, df_train.loc[:, ['1', '2', '3']]], axis=1)
    if 'IHC' in whichFeats:
        df_selecttrain = pd.concat([df_selecttrain, df_train.loc[:, ['HR+HER2+','HR+HER2-','TNBC','HR-HER2+']]], axis=1)
    splits = defineSplits(df_selecttrain, df_ytrain['Status'], random_state)
    print('spllit done')

## MODELS
    ## pCR



    pcr_models = run_all_models(df_selecttrain, yrun, splits, float(rcut),rd=random_state) # 分别是X y splits=cv cut=cut cut是dropcollinear中的corr cut
    pcr_refits = refit_all_models(df_selecttrain, yrun, pcr_models, splits, whichFeats,'c_index',rcut,random_state)
    import os
    stamp = 'submission'
    outdir = 'results_{}_r{}_k{}_rs{}'.format(stamp, rcut, knum,random_state)
    if not os.path.exists('D:/zh/ML/CBCGA_mode2/datav/{}'.format(outdir)):
        os.makedirs('D:/zh/ML/CBCGA_mode2/datav/{}'.format(outdir))
        print('Making: D:/zh/ML/CBCGA_mode2/datav/{}'.format(outdir))
    with open('D:/zh/ML/CBCGA_mode2/datav/{}/{}_pcr_models.pkl'.format(outdir,whichFeats), 'wb') as f:
        pickle.dump(pcr_models, f)
    with open('D:/zh/ML/CBCGA_mode2/datav/{}/{}_pcr_refits.pkl'.format(outdir,whichFeats), 'wb') as f:
        pickle.dump(pcr_refits, f)
    f.close()

if __name__ == "__main__":

    main('RNA_Met_Path_IHC_Clin',0.6,5,111)
