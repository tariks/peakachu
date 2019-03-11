#!/usr/bin/env python
import argparse
from sklearn.externals import joblib
import gc
import pathlib

def main(args):
    import numpy as np
    np.seterr(divide='ignore',invalid='ignore')
    from peakachu import scoreUtils
    pathlib.Path(args.output).mkdir(parents=True, exist_ok=True) 


    import cooler
    Lib = cooler.Cooler(args.path)
    nam = args.path.split('.cool')[0]
    nam = nam.split('/')[-1]
    #resolution = Lib.binsize
    resolution = args.resolution
    chroms = []
    
    model = joblib.load(args.model)
    for key in Lib.chromnames:
        if key.startswith('c'):
            cname=key
        else:
            cname='chr'+key
        X = scoreUtils.Chromosome(Lib.matrix(balance=args.balance, sparse=True).fetch(key).tocsr(),
                                model=model,
                                cname=cname,lower=args.lower,
                                upper=args.upper,res=resolution,
                                width=args.width)
        result,R = X.score()
        X.writeBed(args.output,result,R)

