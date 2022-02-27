#!/usr/bin/env python
import numpy as np
from sklearn.ensemble import RandomForestClassifier as forest
from sklearn.model_selection import GridSearchCV
from sklearn import metrics
from collections import defaultdict, Counter
from scipy import stats
import gc
import random

def buildmatrix(Matrix, coords, width=5, positive=True, stop=5000):
    """
    Generate training set
    :param coords: List of tuples containing coord bins
    :param width: Distance added to center. width=5 makes 11x11 windows
    :return: yields paired positive/negative samples for training
    """
    negcount = 0
    for c in coords:
        x, y = c[0], c[1]
        if y-x < width:
            pass
        else:
            try:
                window = Matrix[x-width:x+width+1,
                                y-width:y+width+1].toarray()
                window[np.isnan(window)] = 0
                if np.count_nonzero(window) < window.size*0.1: # minimum filtering
                    pass
                else:
                    center = window[width, width]
                    p2LL = center/np.mean(window[:width, :width])
                    if positive and p2LL < 0.1:
                        pass
                    else:
                        ranks = stats.rankdata(window, method='ordinal')
                        window = np.r_[window.ravel(), ranks, p2LL]
                        if window.size == 1+2*(2*width+1)**2 and np.all(np.isfinite(window)):
                            if not positive:
                                negcount += 1
                            if negcount >= stop:
                                raise StopIteration
                            yield window
            except:
                pass


def trainRF(X, F, nproc=1):
    """
    :param X: training set from buildmatrix
    :param distances:
    """
    print('input data {} peaks and {} background'.format(
        X.shape[0], F.shape[0]))
    gc.collect()
    params = {}
    params['class_weight'] = ['balanced']
    #params['class_weight'] += [{1: w} for w in range(5,10000,500)]
    params['n_estimators'] = [200]
    params['n_jobs'] = [1]
    params['max_features'] = ['auto']
    params['max_depth'] = [20]
    params['random_state'] = [42]
    #from hellinger_distance_criterion import HellingerDistanceCriterion as hdc
    #h = hdc(1,np.array([2],dtype='int64'))
    params['criterion'] = ['gini', 'entropy']
    #model = forest(**params)
    mcc = metrics.make_scorer(metrics.matthews_corrcoef)
    model = GridSearchCV(forest(), param_grid=params,
                         scoring=mcc, verbose=2, n_jobs=4, cv=5)
    y = np.array([1]*X.shape[0] + [0]*F.shape[0])
    x = np.vstack((X, F))
    model.fit(x, y)
    params = model.best_params_
    print(params)
    print(model.best_score_)
    print('{} peaks {} controls'.format(X.shape[0], F.shape[0]))

    return model.best_estimator_


def parsebed(chiafile, lower=50000, upper=4000000):

    coords = defaultdict(set)
    with open(chiafile) as o:
        for line in o:
            p = line.rstrip().split()

            if 'M' in p[0]:
                continue
            if '_' in p[0]: # remove any unassembled contigs or scaffolds
                continue
            
            s1, e1, s2, e2 = int(p[1]), int(p[2]), int(p[4]), int(p[5])
            if s1 > s2:
                s1, s2 = s2, s1
                e1, e2 = e2, e1
            
            if s2 - s1 > upper:
                continue
            if s2 - s1 < lower:
                continue
            
            # always has prefix "chr", avoid potential bugs
            chrom = 'chr' + p[0].lstrip('chr')
            coords[chrom].add((s1, e1, s2, e2))

    for c in coords:
        coords[c] = sorted(coords[c])

    return coords


def learn_distri_kde(coords, res):

    dis = []
    for c in coords:
        for s1, e1, s2, e2 in coords[c]:
            a = (s1 + e1) // (2 * res)
            b = (s2 + e2) // (2 * res)
            dis.append(b-a)

    lower = min(dis)

    # part 1: same distance distribution as the positive input
    kde = stats.gaussian_kde(dis)

    # part 2: random long-range interactions
    counts, bins = np.histogram(dis, bins=100)
    long_end = int(bins[-1])
    tp = np.where(np.diff(counts) >= 0)[0] + 2
    long_start = int(bins[tp[0]])

    return kde, lower, long_start, long_end


def negative_generating(M, kde, positives, lower, long_start, long_end):

    positives = set(positives)
    N = 3 * len(positives)
    # part 1: kde trained from positive input
    part1 = kde.resample(N).astype(int).ravel()
    part1 = part1[(part1 >= lower) & (part1 <= long_end)]

    # part 2: random long-range interactions
    part2 = []
    pool = np.arange(long_start, long_end+1)
    tmp = np.cumsum(M.shape[0]-pool)
    ref = tmp / tmp[-1]
    for i in range(N):
        r = np.random.random()
        ii = np.searchsorted(ref, r)
        part2.append(pool[ii])

    sample_dis = Counter(list(part1) + part2)

    neg_coords = []
    midx = np.arange(M.shape[0])
    for i in sorted(sample_dis):  # i cannot be zero
        n_d = sample_dis[i]
        R, C = midx[:-i], midx[i:]
        tmp = np.array(M[R, C]).ravel()
        tmp[np.isnan(tmp)] = 0
        mask = tmp > 0
        R, C = R[mask], C[mask]
        pool = set(zip(R, C)) - positives
        sub = random.sample(pool, n_d)
        neg_coords.extend(sub)

    random.shuffle(neg_coords)

    return neg_coords
