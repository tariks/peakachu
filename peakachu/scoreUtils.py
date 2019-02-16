import numpy as np
import os
import gc
from collections import defaultdict
from scipy import sparse
from scipy import stats


class Chromosome():
    def __init__(self,coomatrix,model,lower=1,upper=500,cname='chrm',res=10000,width=5):
        # cLen = coomatrix.shape[0] # seems useless
        R,C = coomatrix.nonzero()
        validmask = np.isfinite(coomatrix.data) & (C-R+1 > lower) & (C-R < upper)
        R,C,data = R[validmask],C[validmask],coomatrix.data[validmask]
        self.M = sparse.csr_matrix((data, (R, C)), shape=coomatrix.shape)
        self.ridx, self.cidx = R, C
        self.chromname = cname
        self.r = res
        self.w=width
        self.model=model


    def getwindow(self,coords):

        """
        Generate training set
        :param Matrix: single chromosome dense array
        :param coords: List of tuples containing coord bins
        :param width: Distance added to center. width=5 makes 11x11 windows
        :return: yields paired positive/negative samples for training
        """
        fts,clist = [],[]
        w2 = int(self.w//2)
        width=self.w
        for c in coords:
            x,y = c[0],c[1]
            distance = abs(y-x)
            try:
                window = self.M[x-width:x+width+1,
                        y-width:y+width+1].toarray()
            except:
                continue
            if np.count_nonzero(window) < window.size*.2:
                pass
            else:
                try:
                    center = window[width,width]
                    ls = window.shape[0]
                    p2LL = center/np.mean(window[ls-1-ls//4:ls,:1+ls//4])
                    indicatar_vars = np.array([p2LL])
                    ranks = stats.rankdata(window,method='ordinal')
                    window = np.hstack((window.flatten(),ranks,indicatar_vars))
                    window = window.flatten()
                    additional=1
     
                    window = window.reshape((1,window.size))
                    if window.size==1+2*(1+2*width)**2 and np.isfinite(window).all():
                        fts.append(window)
                        clist.append(c)
                except:
                    pass
        fts = np.vstack((i for i in fts))
        probas = self.model.predict_proba(fts)[:,1]
        return probas,clist


    def score(self):
        wsize = self.w
        model=self.model
        print('scoring matrix {}'.format(self.chromname))
        print('num candidates {}'.format(self.M.data.size))
        coords = [(r, c) for r, c in zip(self.ridx, self.cidx)]
        p,clist = self.getwindow(coords)
        clist = np.r_[clist]
        pfilter = p >= 0.5
        ri = clist[:,0][pfilter]
        ci = clist[:,1][pfilter]
        result = sparse.csr_matrix((p[pfilter], (ri, ci)), shape=self.M.shape)
        data = np.array(self.M[ri, ci]).ravel()
        self.M = sparse.csr_matrix((data, (ri, ci)), shape=self.M.shape)
        return result, self.M


    def writeBed(self,out,result,R):
        if not os.path.exists(out):
            os.makedirs(out)
        with open(out+'/'+self.chromname+'.bed','w') as o:
            r,c = result.nonzero()
            for i in range(r.size):
                print(self.chromname,r[i]*self.r,(r[i]+1)*self.r,self.chromname,c[i]*self.r,(c[i]+1)*self.r,result[r[i],c[i]],R[r[i],c[i]],sep='\t',file=o)


