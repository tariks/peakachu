#!/usr/bin/env python
# Program to train classifier given a cooler file and
# paired bedfile containing ChIA-PET peaks
# Author: Tarik Salameh

import numpy as np
from random import randint
from sklearn.ensemble import RandomForestClassifier as forest
from sklearn.model_selection import GridSearchCV
from sklearn import metrics
from collections import defaultdict
from scipy import stats
import gc

def buildmatrix(Matrix, coords,width=5,lower=1,positive=True,stop=5000):
    """
    Generate training set
    :param coords: List of tuples containing coord bins
    :param width: Distance added to center. width=5 makes 11x11 windows
    :return: yields paired positive/negative samples for training
    """
    w2 = int(width//2)
    negcount=0
    for c in coords:
        x,y = c[0],c[1]
        if y-x < lower:
            pass
        else:
            try:
                window = Matrix[x-width:x+width+1,
                                y-width:y+width+1].toarray()
                if np.count_nonzero(window) < window.size*.2:
                    pass
                else:
                    center = window[width,width]
                    ls = window.shape[0]
                    p2LL = center/np.mean(window[ls-1-ls//4:ls,:1+ls//4])
                    if positive and p2LL < .1:
                        pass
                    else:
                        indicatar_vars = np.array([p2LL])
                        ranks = stats.rankdata(window,method='ordinal')
                        window = np.hstack((window.flatten(),ranks,indicatar_vars))
                        window = window.flatten()
                        s2 = (1+2*width)**2
                        s2//=2
                        if window.size==1+2*(2*width+1)**2 and np.all(np.isfinite(window)):
                            if not positive:
                                negcount+=1
                            if negcount >=stop:
                                raise StopIteration
                            yield window
            except:
                pass

def trainRF(X,F,nproc=1):
    """
    :param X: training set from buildmatrix
    :param distances:
    """
    print('input data {} peaks and {} background'.format(X.shape[0],F.shape[0]))
    gc.collect()
    params = {}
    #params['class_weight'] = ['balanced',None]
    #params['class_weight'] += [{1: w} for w in range(5,10000,500)]
    params['n_estimators'] = [100]
    params['n_jobs'] = [1]
    params['max_features'] = ['auto']
    params['max_depth'] = [10,20]
    params['random_state'] = [42]
    #from hellinger_distance_criterion import HellingerDistanceCriterion as hdc
    #h = hdc(1,np.array([2],dtype='int64'))
    params['criterion'] = ['entropy','gini']
    #model = forest(**params)
    mcc = metrics.make_scorer(metrics.matthews_corrcoef)
    model = GridSearchCV(forest(),param_grid=params,scoring=mcc,verbose=2,n_jobs=1,cv=5)
    y = np.array([1]*X.shape[0] + [0]*F.shape[0])
    x = np.vstack((X,F))
    model.fit(x,y)
    fts = model.best_estimator_.feature_importances_[:]
    params=model.best_params_
    print(params)
    print(model.best_score_)
    fts = fts.tolist()
    print('{} peaks {} controls'.format(X.shape[0],F.shape[0]))
    return model.best_estimator_


def parsebed(chiafile,res=10000,lower=1):
    coords = defaultdict(list)
    with open(chiafile) as o:
        for line in o:
            s = line.split()
            a,b = float(s[1]),float(s[4])
            a,b = int(a),int(b)
            a//=res
            b//=res
            if b-a > lower and 'Y' not in s[0] and 'X' not in s[0] and 'MT' not in s[0]:
                coords[s[0]].append((a,b))
    return coords

