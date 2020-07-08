import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy import stats

def quantile_norm(X):
    """
    Normalize the columns of X to each have the same distribution.
  
    Given an expression matrix (microarray data, read counts, etc) of M genes
    by N samples, quantile normalization ensures all samples have the same
    spread of data.
  
    The data across each row are averaged to obtain an average column. Each
    column quantile is replaced with the corresponding quantile of the average
    column.
  
    Parameters
    ----------
    X : 2D array of float, shape (M, N)
        The input data, with M rows (genes/features) and N columns (samples).
  
    Returns
    -------
    Xn : 2D array of float, shape (M, N)
        The normalized data.
    """
    # compute the quantiles
    quantiles = np.mean(np.sort(X, axis=0), axis=1)
    # compute the column-wise ranks. Each observation is replaced with its
    # rank in that column: the smallest observation is replaced by 1, the
    # second-smallest by 2, ..., and the largest by M, the number of rows.
    ranks = np.apply_along_axis(stats.rankdata, 0, X)

    # convert ranks to integer indices from 0 to M-1
    rank_indices = ranks.astype(int) - 1

    # index the quantiles for each rank with the ranks matrix
    Xn = quantiles[rank_indices]

    return Xn

def parse_peakachu_loops(fil):

    L = []
    with open(fil, 'r') as source:
        for line in source:
            p = line.rstrip().split()
            key = (p[0], int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5]))
            L.append(key)
    
    return L

def parse_probs(fil):

    L = []
    unique_1 = set()
    unique_2 = set()
    x = []
    y = []
    with open(fil, 'r') as source:
        for line in source:
            p = line.rstrip().split()
            key = (p[0], int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5]))
            if (float(p[6]) > 0) and (float(p[7]) == 0):
                unique_1.add(key)
            elif (float(p[7]) > 0) and (float(p[6]) == 0):
                unique_2.add(key)
            else:
                L.append(key)
                x.append(float(p[6]))
                y.append(float(p[7]))
    x = np.r_[x]
    y = np.r_[y]
    
    return L, unique_1, unique_2, x, y


def work_core(folds, fold1, fold2, union, loop_pool_1, loop_pool_2, thre=0.9):

    # learning the best GMM model
    Pass = False
    for N in range(2, 11):
        model = GaussianMixture(N, covariance_type='full').fit(folds[:,np.newaxis])
        means = model.means_.ravel()
        ri = np.argmax(means)
        probs = model.predict_proba(folds[:,np.newaxis])[:, ri]
        idx = np.where(probs >= thre)[0]
        mask = folds[idx] < 1
        if mask.sum() == 0:
            Pass = True
            break
        else:
            if mask.sum() / mask.size < 0.01:
                Pass = True
                break
      
    unique1 = set()
    unique2 = set()
    if Pass:
        print('Number of Dists: {0}'.format(N))
        probs_1 = model.predict_proba(fold1[:,np.newaxis])[:, ri]
        probs_2 = model.predict_proba(fold2[:,np.newaxis])[:, ri]
        # set1 unique
        idx = np.where(probs_1 >= thre)[0]
        for i in idx:
            key = union[i]
            if (key in loop_pool_1) and (not key in loop_pool_2):
                unique1.add(key)

        # set2 unique
        idx = np.where(probs_2 >= thre)[0]
        for i in idx:
            key = union[i]
            if (key in loop_pool_2) and (not key in loop_pool_1):
                unique2.add(key)
    
    return unique1, unique2, model


loop_pool_1 = set(parse_peakachu_loops(sys.argv[1]))
loop_pool_2 = set(parse_peakachu_loops(sys.argv[2]))
union, unique_1, unique_2, x, y = parse_probs(sys.argv[3])
thre = 0.95
repeat = 50

Xn = quantile_norm(np.r_['1,2,0', x, y])
x, y = Xn.T
fold1 = x / y
fold2 = y / x
folds = np.r_[fold1, fold2]

# to make the model robust to outliers
folds_c = folds.copy()
folds_c.sort()
cumsum = np.cumsum(np.diff(folds_c) < 2)
per_i = np.where(cumsum == np.arange(1, cumsum.size+1))[-0][-1] + 1
per = folds_c[per_i]
folds = folds[folds<=per]


unique1, unique2, model = work_core(folds, fold1, fold2, union, loop_pool_1, loop_pool_2, thre=thre)
for i in range(repeat):
    u1, u2 = work_core(folds, fold1, fold2, union, loop_pool_1, loop_pool_2, thre=thre)[:2]
    unique1 = unique1 & u1
    unique2 = unique2 & u2
    if (not len(unique1)) and (not len(unique2)):
        break

unique1 = unique1 | unique_1
unique2 = unique2 | unique_2
unique1 = sorted(unique1)
unique2 = sorted(unique2)
# output
outfil1 = sys.argv[3].split('.')[0] + '.' + sys.argv[1].replace('.merged.loops', '.unique.loops')
with open(outfil1, 'w') as out:
    for line in unique1:
        tmp = list(map(str, line))
        out.write('\t'.join(tmp)+'\n')

outfil2 = sys.argv[3].split('.')[0] + '.' + sys.argv[2].replace('.merged.loops', '.unique.loops')
with open(outfil2, 'w') as out:
    for line in unique2:
        tmp = list(map(str, line))
        out.write('\t'.join(tmp)+'\n')


# visualize GMM
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
x = np.linspace(0, 5, 500)

logprob = model.score_samples(x[:,np.newaxis])
responsibilities = model.predict_proba(x[:,np.newaxis])
pdf = np.exp(logprob)
pdf_individual = responsibilities * pdf[:,np.newaxis]

ax.hist(folds[folds<5], 50, density=True, histtype='stepfilled', alpha=0.8)
#ax.plot(x, pdf, '-k')
ax.plot(x, pdf_individual, '--k')
ax.set_xlabel('Fold Change', fontsize=18)
ax.set_ylabel('Density', fontsize=18)

#plt.show()
plt.savefig(sys.argv[3].split('.')[0] + '.' + 'Fold-GMM.png', dpi=200, bbox_inches='tight')
plt.close()