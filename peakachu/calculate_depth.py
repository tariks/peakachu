#!/usr/bin/env python
import pathlib

def main(args):
    import argparse, sys
    import numpy as np
    from scipy import sparse
    np.seterr(divide='ignore',invalid='ignore')
    import cooler
    Lib = cooler.Cooler(args.path)
    
    totals=0
    for k in Lib.chromnames[:]:
        intra = np.triu(lib.matrix(balance=False, sparse=False).fetch(k))
        totals += int(intra.sum())
    
    print(totals)


