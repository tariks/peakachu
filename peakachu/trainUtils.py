#!/usr/bin/env python
# Program to train classifier given a cooler file and
# paired bedfile containing ChIA-PET peaks
# Author: Tarik Salameh

import numpy as np
from random import randint
from sklearn.ensemble import RandomForestClassifier as forest
from sklearn.linear_model import SGDClassifier
from sklearn.linear_model import LogisticRegressionCV as LRCV
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import PredefinedSplit
from sklearn.calibration import CalibratedClassifierCV as calib
from sklearn import metrics
from collections import defaultdict
from scipy import stats
import gc

def buildmatrix(Matrix, coords,width=10,lower=1,L=False,chrm='chr1',res=10000,positive=True,stop=2000):
    """
    Generate training set
    :param Matrix: single chromosome dense array
    :param coords: List of tuples containing coord bins
    :param width: Distance added to center. width=5 makes 11x11 windows
    :return: yields paired positive/negative samples for training
    """
    histonedicts = []
    if chrm[0]!='c':
        chrm='chr'+chrm
    w2 = int(width//2)
    negcount=0
    if L:
        import pandas as pd
        for track in L:
            table = pd.read_table(track,usecols=[0,1],index_col=0).loc[chrm]
            table = table.values.T[0].tolist()
            hdict = defaultdict(int)
            for i in table:
                hdict[int(i//res)]+=1
            histonedicts.append(hdict)
    for c in coords:
        x,y = c[0],c[1]
        distance=abs(y-x)
        if y-x < lower:
            pass
        else:
            try:
                # compatible with both numpy.ndarray and CSR/CSC sparse matrix
                window = np.array(Matrix[x-width:x+width+1,
                                         y-width:y+width+1])
                if np.count_nonzero(window) < window.size*.1:
                    pass
                else:
                    center = window[width,width]
                    ls = window.shape[0]
                    p2LL = center/np.mean(window[ls-1-ls//4:ls,:1+ls//4])
                    if positive and p2LL < .5:
                        pass
                    else:
                        indicatar_vars = np.array([p2LL])
                        ranks = stats.rankdata(window,method='ordinal')
                        window = np.hstack((window.flatten(),ranks,indicatar_vars))
                        window = window.flatten()
                        s2 = (1+2*width)**2
                        s2//=2
                        additional=1
                        if len(histonedicts) > 0:
                            for track in histonedicts:
                                window = np.append(window,np.array([track[x],track[y]]))
                                additional+=2
                        #if window.size==s2+additional and np.all(np.isfinite(window)):
                        if window.size==1+2*ranks.size and np.all(np.isfinite(window)):
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
    :return: ExtraTreesClassifier object
    """
    print('input data {} peaks and {} background'.format(X.shape[0],F.shape[0]))
    gc.collect()
    params = {}
    params['class_weight'] = ['balanced',None]
    #params['class_weight'] += [{1: w} for w in range(5,10000,500)]
    params['n_estimators'] = [100]
    params['n_jobs'] = [1]
    params['max_features'] = ['auto']
    params['random_state'] = [42]
    #from hellinger_distance_criterion import HellingerDistanceCriterion as hdc
    #h = hdc(1,np.array([2],dtype='int64'))
    params['criterion'] = ['entropy','gini']
    #model = forest(**params)
    mcc = metrics.make_scorer(metrics.matthews_corrcoef)
    model = GridSearchCV(forest(),param_grid=params,scoring=mcc,verbose=3,n_jobs=1,cv=3)
    y = np.array([1]*X.shape[0] + [0]*F.shape[0])
    x = np.vstack((X,F))
    model.fit(x,y)
    fts = model.best_estimator_.feature_importances_[:]
    params=model.best_params_
    print(params)
    print(model.best_score_)
    fts = fts.tolist()
    print('{} peaks {} controls'.format(X.shape[0],F.shape[0]))
    return model


def parsebed(chiafile,res=10000,lower=1):
    coords = defaultdict(list)
    with open(chiafile) as o:
        for line in o:
            s = line.split()
            a,b = float(s[1]),float(s[4])
            a,b = int(a),int(b)
            a//=res
            b//=res
            if b-a > lower:
                coords[s[0]].append((a,b))
    return coords

