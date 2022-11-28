#!/usr/bin/env python
import numpy as np
from sklearn.ensemble import RandomForestClassifier as forest
from sklearn.model_selection import GridSearchCV
from sklearn import metrics
from collections import defaultdict, Counter
from scipy import stats
import random
from scipy.ndimage import gaussian_filter
from peakachu.utils import image_normalize, distance_normalize, calculate_expected

def buildmatrix(Matrix, coords, w=5):
    """
    Generate training set
    :param coords: List of tuples containing coord bins
    :param width: Distance added to center. width=5 makes 11x11 windows
    :return: yields positive/negative samples for training
    """
    # pre-filter coordinates
    coords = np.r_[coords]
    xi, yi = coords[:,0], coords[:,1]
    mask = (xi - w >= 0) & (yi + w + 1 <= Matrix.shape[0]) & (yi - xi > w)
    xi, yi = xi[mask], yi[mask]
    if xi.size < 10:
        return
    
    # calculate expected values
    maxdis = max([abs(i-j) for i, j in zip(xi, yi)]) + 2*w
    exp_arr = calculate_expected(Matrix, maxdis)
    
    # extract and normalize submatrices surrounding the input coordinates
    seed = np.arange(-w, w+1)
    delta = np.tile(seed, (seed.size, 1))
    xxx = xi.reshape((xi.size, 1, 1)) + delta.T
    yyy = yi.reshape((yi.size, 1, 1)) + delta
    v = np.array(Matrix[xxx.ravel(), yyy.ravel()]).ravel()
    vvv = v.reshape((xi.size, seed.size, seed.size))
    windows, _ = distance_normalize(vvv, exp_arr, xi, yi, w)
    fea = []
    for arr in windows:
        scaled_arr = image_normalize(gaussian_filter(arr, sigma=1, order=0))
        fea.append(scaled_arr.ravel())
    
    return fea

def trainRF(trainset, labels_, nproc=4):
    
    params = {}
    params['class_weight'] = ['balanced', 'balanced_subsample', None]
    params['n_estimators'] = [100]
    params['n_jobs'] = [1]
    params['max_features'] = ['sqrt']
    params['max_depth'] = [15, 20, 25]
    params['criterion'] = ['entropy', 'gini']
    mcc = metrics.make_scorer(metrics.matthews_corrcoef)
    clf = GridSearchCV(forest(), param_grid=params,
                       scoring=mcc, verbose=2, n_jobs=nproc, cv=5)
    clf.fit(trainset, labels_)
    params = clf.best_params_
    print(params)
    print(clf.best_score_)

    return clf.best_estimator_

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
    long_end = max(int(bins[-1]), 5000000//res)
    tp = np.where(np.diff(counts) >= 0)[0] + 2
    long_start = int(bins[tp[0]])

    return kde, lower, long_start, long_end


def negative_generating(M, kde, positives, lower, long_start, long_end):

    positives = set(positives)
    N = 4 * len(positives)
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
        sub = random.sample(pool, min(n_d, len(pool)))
        neg_coords.extend(sub)

    random.shuffle(neg_coords)

    return neg_coords
