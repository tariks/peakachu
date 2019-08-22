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

    if args.path[4:]=='.hic':
        hic=True
    else:
        hic=False
    if not hic:
        import cooler
        Lib = cooler.Cooler(args.path)
        chromosomes=Lib.chromnames[:]
    else:
        import straw
        chromosomes=[str(i) for i in range(1,23)]+['X']
    if cikada in chromosomes:
        ccname=cikada
    else:
        ccname=cikada.split('hr')[1]
    if not hic:
        X = scoreUtils.Chromosome(Lib.matrix(balance=args.balance, sparse=True).fetch(ccname).tocsr(),
                                model=model,
                                cname=cikada,lower=args.lower,
                                upper=args.upper,res=args.resolution,
                                width=args.width)
    elif args.balance:
        X = scoreUtils.Chromosome(straw.straw('KR',args.path,ccname,ccname,'BP',args.resolution),
                                model=model,
                                cname=cikada,lower=args.lower,
                                upper=args.upper,res=args.resolution,
                                width=args.width)
    else:
        X = scoreUtils.Chromosome(straw.straw('NONE',args.path,ccname,ccname,'BP',args.resolution),
                                model=model,
                                cname=cikada,lower=args.lower,
                                upper=args.upper,res=args.resolution,
                                width=args.width)
    result,R = X.score()
    X.writeBed(args.output,result,R)

