#!/usr/env/bin python
import gc
import sys

def main(args):
    import numpy as np
    import pandas as pd
    from collections import defaultdict
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
        idx = (p > args.threshold)
        r,c,p,raw = r[idx],c[idx],p[idx],raw[idx]
        #rawmatrix={(r[i],c[i]): raw[i] for i in range(len(r))}
        matrix={(r[i],c[i]): p[i] for i in range(r.size)}
        D={(r[i],c[i]): raw[i] for i in range(r.size)}
        #P={(r[i],c[i]): p[i] for i in range(r.size)}
        del X
        gc.collect()
        final_list = peakacluster.local_clustering(matrix,res=res)
        final_list = [i[0] for i in final_list]
        r = [i[0] for i in final_list]
        c = [i[1] for i in final_list]
        p = np.array([matrix.get((r[i],c[i])) for i in range(len(r))])
        for i in range(len(r)):
            P = matrix.get((r[i],c[i]))
            print(chrom,r[i]*res,r[i]*res+res,chrom,
                  c[i]*res,c[i]*res+res,
                  P,sep='\t')
