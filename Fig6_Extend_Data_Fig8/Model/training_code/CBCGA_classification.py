from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest

from sklearn.pipeline import Pipeline
from sklearn.model_selection import RandomizedSearchCV,GridSearchCV

from sklearn.base import BaseEstimator, TransformerMixin, RegressorMixin, clone

from scipy.stats import spearmanr

from copy import deepcopy
import numpy as np
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored

from sksurv.preprocessing import OneHotEncoder
from sksurv.ensemble import RandomSurvivalForest
from sksurv.svm import FastKernelSurvivalSVM, FastSurvivalSVM





def optimise_cox_featsel(X,y,cut,cv=5):
    cox = CoxnetSurvivalAnalysis()
    pipe = Pipeline([('model', cox)])
    param_grid = {'model__l1_ratio':np.arange(0.1,1.1,0.1),}
    search = GridSearchCV(pipe, param_grid, cv=cv, scoring=score_survival_model,n_jobs=50)
    search.fit(X,y)
    print(search.best_params_)
    return search

def optimise_rf_featsel(X,y,cut,cv=5,rd=None):
    rf=RandomSurvivalForest(random_state=rd)
    pipe = Pipeline([('rf', rf)])
    param_grid = {"rf__max_depth": [3, None],
                  "rf__n_estimators": [5, 10, 25, 50, 100],
                  "rf__max_features": [0.05, 0.1, 0.2, 0.5, 0.7],
                  "rf__min_samples_split": [2, 3, 6, 10, 12, 15]
                  }
    search = GridSearchCV(pipe, param_grid, cv=cv, scoring=score_survival_model,n_jobs=50)
    search.fit(X,y)
    print(search.best_params_)
    return search

def optimise_svc_featsel(X,y,cut,cv=5,rd=None):
    svc=FastSurvivalSVM(random_state=rd)
    pipe = Pipeline([('svc', svc)])
    param_grid = { 'svc__alpha': 2. ** np.arange(-12, 13, 2)}
    search = GridSearchCV(pipe, param_grid, cv=cv, scoring=score_survival_model, n_jobs=50)
    search.fit(X, y)
    print(search.best_params_)
    return search

def plot_and_refit(X, y, model, cv, label='c_index', prefix='someresponse',feats='features',rcut=None,random_state=None):
    ypreds = []
    yreals = []
    ypreds_cv = []
    yreals_cv = []
    cv_models = []
    c_index = []
    for i,(tr,ts) in enumerate(cv):
        model.fit(X.iloc[tr,:], y[tr])
        cv_models.append(deepcopy(model))
        y_pred = model.predict(X.iloc[ts,:])
        ytest = y[ts]

        # Precision
        ypreds.extend(y_pred)
        yreals.extend(ytest)
        ypreds_cv.append(y_pred)
        yreals_cv.append(ytest)
        c = concordance_index_censored(ytest['Status'],ytest['Survival_in_days'],y_pred)[0]
        c_index.append(c)
        #c_index
    mean_c_index=np.mean(c_index)
    f = open('D:/zh/ML/CBCGA_mode2/output/'+'output_all_rcut{}_rd{}.csv'.format(rcut,random_state), 'a')
    f.write('{},cv,{},{},{},{},{}\n'.format(feats,prefix,label,mean_c_index,rcut,random_state))
    print('{},cv,{},{},{},{},{}\n'.format(feats,prefix,label,mean_c_index,rcut,random_state))
    f.close()

    model.fit(X,y)
    return [model, cv_models]

def run_all_models(X,y,splits,cut,rd=None):
    cox_result = optimise_cox_featsel(X,y,cut=cut,cv=splits)
    svc_result = optimise_svc_featsel(X,y,cut=cut,cv=splits,rd=rd)
    rf_result = optimise_rf_featsel(X,y,cut=cut,cv=splits,rd=rd)
    averaged_models = AveragingModels(models = (cox_result,svc_result,rf_result)) #svc使用predict函数就可以了
    results = {}
    results['lr'] = cox_result
    results['svc'] = svc_result
    results['rf'] = rf_result
    results['avg'] = averaged_models
    return results

def refit_all_models(X,y,results,splits,whichFeats,criterion,rcut,random_state):
    refit = {}
    for model in results.keys():
        try:
            refit[model] = plot_and_refit(X,y,results[model].best_estimator_,splits,label=model,prefix=criterion,feats=whichFeats,rcut=rcut,random_state=random_state)
        except:
            refit[model] = plot_and_refit(X,y,results[model],splits,label=model,prefix=criterion,feats=whichFeats,rcut=rcut,random_state=random_state)
    return refit



class SelectAtMostKBest(SelectKBest):
    def _check_params(self, X, y):
        if not (self.k == "all" or 0 <= self.k <= X.shape[1]):
            # set k to "all" (skip feature selection), if less than k features are available
            self.k = "all"

def score_survival_model(model, X, y):
    prediction = model.predict(X)
    result = concordance_index_censored(y['Status'], y['Survival_in_days'], prediction)
    return result[0]

def fit_and_score_features(X, y):
    n_features = X.shape[1]
    scores = np.empty(n_features)
    m = CoxPHSurvivalAnalysis()
    for j in range(n_features):
        Xj = X[:, j:j+1]
        m.fit(Xj, y)
        scores[j] = m.score(Xj, y)
    return scores

class AveragingModels(BaseEstimator, RegressorMixin, TransformerMixin):
    def __init__(self, models):
        self.models = models

    # we define clones of the original models to fit the data in
    def fit(self, X, y):
        self.models_ = [clone(x.best_estimator_) for x in self.models]

        # Train cloned base models
        for model in self.models_:
            model.fit(X, y)

        return self

    #Now we do the predictions for cloned models and average them
    def predict(self, X):
        from sklearn import preprocessing
        zscore = preprocessing.StandardScaler()
        prediction_cox=self.models[0].predict(X)
        prediction_svc=self.models[1].predict(X)
        prediction_rf=self.models[2].predict(X)
        prediction_cox_scaler=zscore.fit_transform(prediction_cox.reshape(-1,1))
        prediction_svc_scaler=zscore.fit_transform(prediction_svc.reshape(-1,1))
        prediction_rf_scaler=zscore.fit_transform(prediction_rf.reshape(-1,1))
        predict_avg=np.mean(np.concatenate((prediction_cox_scaler,prediction_svc_scaler,prediction_rf_scaler),axis=1),axis=1)

        return predict_avg

