#!/usr/bin/env python
import argparse, pickle
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
    resolution = 10000
    chroms = []

    mod = args.model
    with open(mod,'rb') as o:
        model = pickle.load(o)
    for key in Lib.chromnames:
        if key.startswitch('c'):
            cname=key
        else:
            cname='chr'+key
        X = scoreUtils.Chromosome(Lib.matrix(balance=True, sparse=True).fetch(key).tocsr(),
                                model=model,
                                cname=cname,lower=args.lower,
                                upper=args.upper,res=resolution,
                                width=args.width)
        result,R = X.score()
        X.writeBed(args.output,result,R)

