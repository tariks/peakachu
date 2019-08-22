#!/usr/bin/env python
import pathlib

def main(args):
    import argparse, sys
    import numpy as np
    from scipy import sparse
    np.seterr(divide='ignore',invalid='ignore')
    if args.path[-3:]=='.hic':
        hic=True
    else:
        hic=False

    totals=0
    if not hic:
        import cooler
        Lib = cooler.Cooler(args.path)
    
        for k in Lib.chromnames[:]:
            intra = np.triu(Lib.matrix(balance=False, sparse=False).fetch(k))
            totals += int(intra.sum())
    else:
        import straw
        chromosomes = list(range(1,23))+['X'] 
        for k in chromosomes:
            intra = straw.straw('NONE',args.path,k,k,'BP',10000)
            totals += int(intra.sum())

    print(totals)


