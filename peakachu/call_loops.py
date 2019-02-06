#!/usr/env/bin python
import gc

def main(args):
    import numpy as np
    import pandas as pd
    from scipy import ndimage
    from peakachu import peakacluster
    res = 10
    x = pd.read_table(args.infile,index_col=0,
                      usecols=[0,1,4,6,7],header=None)
    chromosomes = list(set(x.index))

    for chrom in chromosomes:
        X = x.loc[chrom].values
        r = X[:,0].astype(int)//res
        c = X[:,1].astype(int)//res
        p = X[:,2].astype(float)
        m = (p>.6)
        raw = X[:,3].astype(float)
        r,c,p,raw = r[m],c[m],p[m],raw[m]
        maxy = np.max(c)+1
        matrix = np.zeros((maxy,maxy))
        rawmatrix = np.zeros((maxy,maxy))
        matrix[r,c]+=p
        rawmatrix[r,c]+=raw
        screen = np.zeros((21,21))
        screen+=np.tril(np.ones((21,21)),2)
        screen-=np.tril(screen,-8)
        matrix+=np.triu(matrix.T,-9)
        percentile_filter = ndimage.percentile_filter(matrix,50,footprint=screen) 
        percentile_filter-=np.tril(percentile_filter,5)
        matrix[matrix<percentile_filter]=0
        del percentile_filter
        gc.collect()
        r,c = matrix.nonzero()
        p = matrix[r,c]
        del X
        gc.collect()
        idx = np.argsort(p)
        idx = idx[-20000:]
        r,c = r[idx],c[idx]
        rawmatrix[matrix==0]=0
        raw = rawmatrix[r,c]
        idx = np.argsort(raw)
        idx = idx[-15000:]
        D = {(r[i],c[i]): raw[i] for i in idx}
        gc.collect()
        final_list = peakacluster.local_clustering(D,res)
        final_list = [i[0] for i in final_list]
        r = [i[0] for i in final_list]
        c = [i[1] for i in final_list]
        if len(r) > 2000:
            p=matrix[r,c]
            sorted_index=np.argsort(p)
            r = [r[i] for i in sorted_index[-2000:]]
            c = [c[i] for i in sorted_index[-2000:]]
        for i in range(len(r)):
            P=matrix[r[i],c[i]]
            print(chrom,r[i]*res,r[i]*res+res,chrom,
                  c[i]*res,c[i]*res+res,
                  P,rawmatrix[r[i],c[i]],
                  sep='\t')
