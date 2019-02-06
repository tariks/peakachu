#!/usr/bin/env python
import argparse, pickle
import gc
import pathlib

def main(args):
    import numpy as np
    np.seterr(divide='ignore',invalid='ignore')
    from peakachu import scoreUtils
    pathlib.Path(args.output).mkdir(parents=True, exist_ok=True) 
    cikada=args.model
    cikada=cikada.split('.pk')[0]
    cikada='ch'+cikada.split('ch')[-1]


    mod = args.model
    with open(mod,'rb') as o:
        model = pickle.load(o)
    import cooler
    Lib = cooler.Cooler(args.path)
    if cikada in Lib.chromnames:
        ccname=cikada
    else:
        ccname=cikada.split('hr')[1]
    #resolution = Lib.binsize
    resolution = 10000
    chroms = []

    X = scoreUtils.Chromosome(Lib.matrix(balance=True, sparse=False).fetch(ccname),
                                model=model,
                                cname=cikada,lower=args.lower,
                                upper=args.upper,res=resolution,
                                width=args.width)
    result,R = X.score()
    X.writeBed(args.output,result,R)

