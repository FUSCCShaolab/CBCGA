from sklearn.feature_selection import SelectKBest
from sklearn.base import BaseEstimator, TransformerMixin, RegressorMixin, clone
import numpy as np
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.metrics import concordance_index_censored

def average(data):
    return sum(data) / len(data)


def bootstrap_cindex(model, X_test, y_test, nsamples=1000):
    cindex_values = []
    for b in range(nsamples):
        idx = np.random.randint(X_test.shape[0], size=X_test.shape[0])
        pred=model.predict(X_test.iloc[idx])
        cindex=concordance_index_censored(y_test[idx]['Status'],y_test[idx]['Survival_in_days'],pred)[0]
        cindex_values.append(cindex)
    return np.percentile(cindex_values, (2.5, 97.5))



def final_test(X, y, model, label='Response', prefix='someresponse',rcut=None,randomseed=None,whichfeat=None,knum=None):

    y_pred = model.predict(X)
    ytest_cindex=concordance_index_censored(y['Status'],y['Survival_in_days'],y_pred)[0]
    cindex_low95,cindex_up95=bootstrap_cindex(model,X,y,nsamples=1000)
    f = open('D:/zh/ML/CBCGA_mode2/testoutput/r{}_knum{}_rs{}.csv'.format(rcut,knum,randomseed), 'a')
    print(whichfeat,'_',ytest_cindex,'_',cindex_low95,'_',cindex_up95)
    f.write('test,{},{},{},{},{}\n'.format(whichfeat,label,ytest_cindex,cindex_low95,cindex_up95))
    f.close()



def test_all_models(X,y,results,criterion,rcut=None,randomseed=None,whichfeat=None,knum=None):
    test_result = {}
    mode=[]
    for model in results.keys():
        mode.append(model)
    test_result[mode[2]] = final_test(X,y,results[mode[2]],label=model,prefix=criterion,rcut=rcut,randomseed=randomseed,whichfeat=whichfeat,knum=knum)
    return test_result

class AveragingModels(BaseEstimator, RegressorMixin, TransformerMixin):
    def __init__(self, models):
        self.models = models

    def fit(self, X, y):
        self.models_ = [clone(x.best_estimator_) for x in self.models]

        # Train cloned base models
        for model in self.models_:
            model.fit(X, y)

        return self

    def predict(self, X):
        from sklearn import preprocessing
        zscore = preprocessing.StandardScaler()
        prediction_cox=self.models[0].predict(X)
        prediction_svc=self.models[1].predict(X)
        prediction_rf=self.models[2].predict(X)
        prediction_cox_scaler=zscore.fit_transform(prediction_cox.shape(-1,1))
        prediction_svc_scaler=zscore.fit_transform(prediction_svc.shape(-1,1))
        prediction_rf_scaler=zscore.fit_transform(prediction_rf.shape(-1,1))
        predict_avg=np.mean(np.concatenate((prediction_cox_scaler,prediction_svc_scaler,prediction_rf_scaler),axis=1),axis=1)

        return predict_avg



### Custom class inspired by:
### https://stackoverflow.com/questions/25250654/how-can-i-use-a-custom-feature-selection-function-in-scikit-learns-pipeline
### https://ramhiser.com/post/2018-04-16-building-scikit-learn-pipeline-with-pandas-dataframe/
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

### Inspired by: https://stackoverflow.com/questions/29412348/selectkbest-based-on-estimated-amount-of-features/29412485
class SelectAtMostKBest(SelectKBest):
    def _check_params(self, X, y):
        if not (self.k == "all" or 0 <= self.k <= X.shape[1]):
            # set k to "all" (skip feature selection), if less than k features are available
            self.k = "all"

def fit_and_score_features(X, y):
    n_features = X.shape[1]
    scores = np.empty(n_features)
    m = CoxPHSurvivalAnalysis()
    for j in range(n_features):
        Xj = X[:, j:j+1]
        m.fit(Xj, y)
        scores[j] = m.score(Xj, y)
    return scores