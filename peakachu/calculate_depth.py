#!/usr/bin/env python
import pathlib

def main(args):
    import argparse, sys
    import numpy as np
    from scipy import sparse
    np.seterr(divide='ignore',invalid='ignore')
    import cooler
    Lib = cooler.Cooler(args.path)
    #resolution = Lib.binsize
    resolution = 10000
    sums=0
    totals=0
    for k in Lib.chromnames[:]:
        x=Lib.matrix(balance=False,sparse=False).fetch(k).tocsr()
        sums+=sparse.triu(x,k=2).sum()
        sums-=sparse.triu(x,k=300).sum()
        totals+=sparse.triu(x).sum()

    print(sums)
    print(totals)


