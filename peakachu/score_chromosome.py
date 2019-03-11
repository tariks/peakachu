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
    cikada=args.model
    cikada=cikada.split('.pk')[0]
    cikada='ch'+cikada.split('ch')[-1]


    model = joblib.load(args.model)

    import cooler
    Lib = cooler.Cooler(args.path)
    if cikada in Lib.chromnames:
        ccname=cikada
    else:
        ccname=cikada.split('hr')[1]
    #resolution = Lib.binsize
    resolution = args.resolution

    X = scoreUtils.Chromosome(Lib.matrix(balance=args.balance, sparse=True).fetch(ccname).tocsr(),
                                model=model,
                                cname=cikada,lower=args.lower,
                                upper=args.upper,res=resolution,
                                width=args.width)
    result,R = X.score()
    X.writeBed(args.output,result,R)

