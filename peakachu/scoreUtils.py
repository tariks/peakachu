from collections import defaultdict
import numpy as np
from scipy import sparse
from scipy import stats
from scipy.ndimage import gaussian_filter
from peakachu.utils import image_normalize, distance_normalize, calculate_expected
#from statsmodels.sandbox.stats.multicomp import multipletests

class Chromosome():
    def __init__(self, M, model, raw_M=None, weights=None,
        lower=6, upper=300, cname='chrm', res=10000, width=5):
        
        lower = max(lower, width+1)
        upper = min(upper, M.shape[0]-2*width)
        # calculate expected values
        if weights is None:
            self.exp_arr = calculate_expected(M, upper+2*width, raw=True)
            if M is raw_M:
                self.background = self.exp_arr
            else:
                self.background = calculate_expected(raw_M, upper+2*width, raw=True)
        else:
            self.exp_arr = calculate_expected(M, upper+2*width, raw=False)
            self.background = self.exp_arr
        
        self.raw_M = raw_M
        self.weights = weights

        # lower down the memory usage
        R, C = M.nonzero()
        validmask = np.isfinite(M.data) & (C-R > (-2*width)) & (C-R < (upper+2*width))
        R, C, data = R[validmask], C[validmask], M.data[validmask]
        self.M = sparse.csr_matrix((data, (R, C)), shape=M.shape)
        self.get_candidate(lower, upper)
        self.chromname = cname
        self.r = res
        self.w = width
        self.model = model
    
    def get_candidate(self, lower, upper):

        x_arr = np.array([], dtype=int)
        y_arr = np.array([], dtype=int)
        p_arr = np.array([], dtype=float)
        idx = np.arange(self.raw_M.shape[0])
        for i in range(lower, upper+1):
            diag = self.raw_M.diagonal(i)
            e = self.background[i]
            if (diag.size > 0) and (e > 0):
                xi = idx[:-i]
                yi = idx[i:]
                if self.weights is None:
                    exp = np.ones(diag.size, dtype=float) * e
                else:
                    b1 = self.weights[:-i]
                    b2 = self.weights[i:]
                    exp = np.ones(diag.size, dtype=float) * e / (b1 * b2)
                
                Poiss = stats.poisson(exp)
                pvalues = Poiss.sf(diag)
                mask = (diag > 0) & np.isfinite(pvalues)
                x_arr = np.r_[x_arr, xi[mask]]
                y_arr = np.r_[y_arr, yi[mask]]
                p_arr = np.r_[p_arr, pvalues[mask]]
        
        #qvalues = multipletests(p_arr, method = 'fdr_bh')[1]
        mask = p_arr < 0.01
        self.ridx, self.cidx = x_arr[mask], y_arr[mask]

    def getwindow(self, coords):
        
        w = self.w
        coords = np.r_[coords]
        xi, yi = coords[:,0], coords[:,1]
        mask = (xi - w >= 0) & (yi + w + 1 <= self.M.shape[0])
        xi, yi = xi[mask], yi[mask]
        seed = np.arange(-w, w+1)
        delta = np.tile(seed, (seed.size, 1))
        xxx = xi.reshape((xi.size, 1, 1)) + delta.T
        yyy = yi.reshape((yi.size, 1, 1)) + delta
        v = np.array(self.M[xxx.ravel(), yyy.ravel()]).ravel()
        vvv = v.reshape((xi.size, seed.size, seed.size))
        windows, clist = distance_normalize(vvv, self.exp_arr, xi, yi, w)
        fea = []
        for arr in windows:
            tmp = gaussian_filter(arr, sigma=1, order=0)
            scaled_arr = image_normalize(tmp)
            fea.append(scaled_arr.ravel())

        fea = np.r_[fea]
        clist = np.r_[clist]

        return fea, clist

    def score(self, thre=0.5):

        print('scoring matrix {}'.format(self.chromname))
        print('number of candidates {}'.format(self.ridx.size))
        total_coords = [(r, c) for r, c in zip(self.ridx, self.cidx)]
        ri = np.r_[[]]
        ci = np.r_[[]]
        prob_pool = np.r_[[]]
        # to lower down the memory usage
        batch_size = 100000
        for t in range(0, len(total_coords), batch_size):
            coords = total_coords[t:t+batch_size]
            fea, clist = self.getwindow(coords)
            if fea.shape[0] > 1:
                p = self.model.predict_proba(fea)[:, 1]
                pfilter = p > thre
                ri = np.r_[ri, clist[:, 0][pfilter]]
                ci = np.r_[ci, clist[:, 1][pfilter]]
                prob_pool = np.r_[prob_pool, p[pfilter]]
        ri = ri.astype(int)
        ci = ci.astype(int)
        
        
        result = sparse.csr_matrix((prob_pool, (ri, ci)), shape=self.M.shape)
        if ri.size > 0:
            data = np.array(self.M[ri, ci]).ravel()
            self.M = sparse.csr_matrix((data, (ri, ci)), shape=self.M.shape)
        else:
            self.M = result

        return result, self.M

    def writeBed(self, outfil, prob_csr, raw_csr):
        
        with open(outfil, 'a') as out:
            r, c = prob_csr.nonzero()
            for i in range(r.size):
                line = [self.chromname, r[i]*self.r, (r[i]+1)*self.r,
                        self.chromname, c[i]*self.r, (c[i]+1)*self.r,
                        prob_csr[r[i],c[i]], raw_csr[r[i],c[i]]]
                out.write('\t'.join(list(map(str, line)))+'\n')
