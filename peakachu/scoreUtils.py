from collections import defaultdict
import numpy as np
from scipy import sparse
from scipy import stats

class Chromosome():
    def __init__(self, coomatrix, model, lower=6, upper=300, cname='chrm', res=10000, width=5):
        
        lower = max(lower, width+1)
        R, C = coomatrix.nonzero()
        validmask = np.isfinite(coomatrix.data) & (C-R > (-2*width)) & (C-R < (upper+2*width))
        R, C, data = R[validmask], C[validmask], coomatrix.data[validmask]
        self.M = sparse.csr_matrix((data, (R, C)), shape=coomatrix.shape)
        self.get_candidate(lower, upper)
        self.chromname = cname
        self.r = res
        self.w = width
        self.model = model
    
    def get_candidate(self, lower, upper):

        x_arr = np.array([], dtype=int)
        y_arr = np.array([], dtype=int)
        idx = np.arange(self.M.shape[0])
        for i in range(lower, upper+1):
            xi = idx[:-i]
            yi = idx[i:]
            diag = self.M.diagonal(i)
            if diag.size > 0:
                e = diag.mean()
                if e > 0:
                    diag = diag / e
                    mask = diag > 1
                    x_arr = np.r_[x_arr, xi[mask]]
                    y_arr = np.r_[y_arr, yi[mask]]
        
        self.ridx, self.cidx = x_arr, y_arr

    def getwindow(self, coords):
        
        fea, clist = [], []
        w = self.w
        for x, y in coords:
            if (x - w < 0) or (y + w + 1 > self.M.shape[0]):
                continue

            window = self.M[x-w:x+w+1, y-w:y+w+1].toarray()
            if np.count_nonzero(window) < window.size*0.1:
                continue
            
            if np.mean(window[:w, :w]) > 0:
                ranks = stats.rankdata(window, method='ordinal')
                center = window[w, w]
                p2LL = center/np.mean(window[:w, :w])
                window = np.r_[window.ravel(), ranks, p2LL]
                if window.size == 1+2*(1+2*w)**2:
                    fea.append(window)
                    clist.append((x, y))

        fea = np.r_[fea]
        clist = np.r_[clist]

        return fea, clist

    def score(self, thre=0.5):

        print('scoring matrix {}'.format(self.chromname))
        print('number of candidates {}'.format(self.ridx.size))
        coords = [(r, c) for r, c in zip(self.ridx, self.cidx)]
        fea, clist = self.getwindow(coords)
        p = self.model.predict_proba(fea)[:, 1]
        pfilter = p > thre
        ri = clist[:, 0][pfilter]
        ci = clist[:, 1][pfilter]
        result = sparse.csr_matrix((p[pfilter], (ri, ci)), shape=self.M.shape)
        data = np.array(self.M[ri, ci]).ravel()
        self.M = sparse.csr_matrix((data, (ri, ci)), shape=self.M.shape)

        return result, self.M

    def writeBed(self, outfil, prob_csr, raw_csr):
        
        with open(outfil, 'a') as out:
            r, c = prob_csr.nonzero()
            for i in range(r.size):
                line = [self.chromname, r[i]*self.r, (r[i]+1)*self.r,
                        self.chromname, c[i]*self.r, (c[i]+1)*self.r,
                        prob_csr[r[i],c[i]], raw_csr[r[i],c[i]]]
                out.write('\t'.join(list(map(str, line)))+'\n')
