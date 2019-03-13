#!/usr/env/bin python
import gc
import sys

def main(args):
    import numpy as np
    import pandas as pd
    from scipy import ndimage
    from peakachu import peakacluster
    res = args.resolution
    x = pd.read_table(args.infile,index_col=0,
                      usecols=[0,1,4,6,7],header=None)
    chromosomes = list(set(x.index))

    for chrom in chromosomes:
        X = x.loc[chrom].values
        r = X[:,0].astype(int)//res
        c = X[:,1].astype(int)//res
        p = X[:,2].astype(float)
        raw = X[:,3].astype(float)
        d = c-r
        unique_d=list(set(d.tolist()))
        idx = (p > args.threshold)
        r,c,p,raw,d = r[idx],c[idx],p[idx],raw[idx],d[idx]
        #rawmatrix={(r[i],c[i]): raw[i] for i in range(len(r))}
        matrix={(r[i],c[i]): p[i] for i in range(len(r))}
        count=p.size
        while count > 10000:
            D={}
            for distance in unique_d:
                dx=(d==distance)
                dr,dc,dp,draw=r[dx],c[dx],p[dx],raw[dx]
                dx=(dp>np.median(dp))
                dr,dc,dp,draw=dr[dx],dc[dx],dp[dx],draw[dx]
                for i in range(dr.size):
                    D[(dr[i],dc[i])]=draw[i]
            count=len(D.keys())

        del X
        gc.collect()
        #idx = np.argsort(p)
        #idx = idx[0-idx.size//2:]
        #idx = np.argsort(raw[idx])
        #idx = idx[0-idx.size//2:]
        #D = {(r[i],c[i]): rawmatrix.get((r[i],c[i])) for i in idx}
        final_list = peakacluster.local_clustering(D,res=res)
        final_list = [i[0] for i in final_list]
        r = [i[0] for i in final_list]
        c = [i[1] for i in final_list]
        p = np.array([matrix.get((r[i],c[i])) for i in range(len(r))])
        if len(r) > 2000:
            sorted_index=np.argsort(p)
            r = [r[i] for i in sorted_index[-2000:]]
            c = [c[i] for i in sorted_index[-2000:]]
        for i in range(len(r)):
            P = matrix.get((r[i],c[i]))
            print(chrom,r[i]*res,r[i]*res+res,chrom,
                  c[i]*res,c[i]*res+res,
                  P,sep='\t')
