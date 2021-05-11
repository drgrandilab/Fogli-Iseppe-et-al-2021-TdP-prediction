## This code runs a Recursive Feature Elimination algorithm with hyperparameters optimization.
## At each step, a Support Vector Machine is trained, and the regularization parameter C (1/lambda) is optimized using the hyperopt package (bayesian optimizer).
## The cost function to minimize is the Matthews' Correlation Coefficient calculated on Leave-One-Out cross-validation changed of sign.



## import packages

import numpy as np
import pandas as pd
import os
from joblib import Parallel, delayed
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import f1_score, matthews_corrcoef, roc_auc_score
from hyperopt import hp, atpe, fmin, STATUS_OK, Trials

## define custom functions

def create_feature_ranking(feat_names, feat_coef):

	## sort features by their coefficients

    df = sorted([(abs(feat_coef[i]), feat_names[i]) for i in range(len(feat_names))], reverse=True)
    return(df)


def Regularized_svm_rfe(X_stand, y, C):

	## calculate MCC through LOO with a given value of C

    svm = SVC(C=C, kernel='linear', max_iter=10e5, probability=True, random_state=42)
    y_pred = []
    loo = LeaveOneOut()

    for train_index, test_index in loo.split(X_stand):
        X_train, X_test = X_stand[train_index,:], X_stand[test_index,:]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        svm.fit(X_train, y_train)
        y_pred.append(svm.predict(X_test)[0])

    return({'loss': -matthews_corrcoef(y, y_pred), 'status': STATUS_OK})


def run_optimization(X, y, n_iter, results):

	## identify best performing value of C with bayesian optimization and return the updated dataset with the least important feature removed
    
    X_stand = StandardScaler().fit_transform(X)
    
    trials = Trials()
        
    best = fmin(fn = lambda C: Regularized_svm_rfe(X_stand, y, C),
                space = hp.uniform('C', 0.001, 100), 
                algo=atpe.suggest, 
                max_evals = n_iter,
                trials=trials)
    
    svm = SVC(C=best['C'], 
              kernel='linear', 
              max_iter=10e5,
              random_state=42)
    
    svm.fit(X_stand, y)
    
    feat_rank = create_feature_ranking(X.columns, svm.coef_[0])
    
    results.append([int(len(X.columns)),
                    trials.best_trial['misc']['vals']['C'][0],
                    -trials.best_trial['result']['loss'],
                    list(X.columns)])
    
    feat_selection = [x for (i,x) in feat_rank[:-1]]
    
    return(X[feat_selection])


def repeat(data, file, trial):

	# run the entire RFE algorithm for 'n_iter' times and save the results in a csv file

    n_iter = 200

    X = data.iloc[:, 2:-1].copy()
    y = data['TdP'].copy()
    results = []
    max_feats = X.shape[1]

    for i in range(max_feats):

        print(max_feats-i)
        
        X = run_optimization(X,y,n_iter,results)
        
    output = pd.DataFrame(results, columns=['n_features', 'C', 'MCC', 'features'])

    output.to_csv('./results/%s_svm_optimization_%i_iter_trial_%i.csv' % (file[4:-13], n_iter, trial))


def run_rfe_with_bayes(file):

	# read the initial dataset and parallelize the code

    data = pd.read_csv('./data/'+file, index_col=0)

    Parallel(n_jobs=10)(delayed(repeat)(data, file ,trial) for trial in range(1)) ### n_jobs is the number of cores that you want to use for the parallelization



##############################


file_list = [i for i in os.listdir('./data/')]

for i in file_list:
    run_rfe_with_bayes(i)
