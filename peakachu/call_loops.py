#!/usr/env/bin python


def main(args):

    import gc
    import sys
    import numpy as np
    import pandas as pd
    from collections import defaultdict
    from peakachu import peakacluster

    res = args.resolution
    x = {}
    with open(args.infile, 'r') as source:
        for line in source:
            p = line.rstrip().split()
            chrom = p[0]
            if float(p[6]) > args.threshold:
                if not chrom in x:
                    x[chrom] = []
                x[chrom].append([int(p[1]), int(p[4]), float(p[6]), float(p[7])])
    for c in x:
        x[c] = np.r_[x[c]]

    for chrom in x:
        X = x[chrom]
        r = X[:, 0].astype(int)//res
        c = X[:, 1].astype(int)//res
        p = X[:, 2].astype(float)
        raw = X[:, 3].astype(float)
        d = c-r
        tmpr, tmpc, tmpp, tmpraw, tmpd = r, c, p, raw, d
        #rawmatrix={(r[i],c[i]): raw[i] for i in range(len(r))}
        matrix = {(r[i], c[i]): p[i] for i in range(len(r))}
        count = 40001
        while count > 40000:
            D = defaultdict(float)
            P = defaultdict(float)
            unique_d = list(set(tmpd.tolist()))
            for distance in unique_d:
                dx = (tmpd == distance)
                dr, dc, dp, draw = tmpr[dx], tmpc[dx], tmpp[dx], tmpraw[dx]
                dx = (dp > np.percentile(dp, 10))
                dr, dc, dp, draw = dr[dx], dc[dx], dp[dx], draw[dx]
                for i in range(dr.size):
                    D[(dr[i], dc[i])] += draw[i]
                    P[(dr[i], dc[i])] += dp[i]
            count = len(D.keys())
            tmpr = np.array([i[0] for i in P.keys()])
            tmpc = np.array([i[1] for i in P.keys()])
            tmpp = np.array([P.get(i) for i in P.keys()])
            tmpraw = np.array([D.get(i) for i in P.keys()])
            tmpd = tmpc-tmpr

        del X
        gc.collect()
        final_list = peakacluster.local_clustering(D, res=res)
        final_list = [i[0] for i in final_list]
        r = [i[0] for i in final_list]
        c = [i[1] for i in final_list]
        p = np.array([matrix.get((r[i], c[i])) for i in range(len(r))])
        if len(r) > 7000:
            sorted_index = np.argsort(p)
            r = [r[i] for i in sorted_index[-7000:]]
            c = [c[i] for i in sorted_index[-7000:]]
        for i in range(len(r)):
            P = matrix.get((r[i], c[i]))
            print(chrom, r[i]*res, r[i]*res+res, chrom,
                  c[i]*res, c[i]*res+res,
                  P, sep='\t')
