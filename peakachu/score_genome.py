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

    if args.path[-4:] == '.hic':
        hic=True
        import cooler
        Lib = cooler.Cooler(args.path)
        chromosomes=Lib.chromnames[:]
        nam = args.path.split('.cool')[0]
    else:
        hic=False
        import straw
        chromosomes=[str(i) for i in range(1,23)]+['X']
        nam = args.path.split('.hic')[0]
    nam = nam.split('/')[-1]
    model = joblib.load(args.model)
    for key in chromosomes:
        if key.startswith('c'):
            cname=key
        else:
            cname='chr'+key
        if not hic:
            X = scoreUtils.Chromosome(Lib.matrix(balance=args.balance, sparse=True).fetch(key).tocsr(),
                                model=model,
                                cname=cname,lower=args.lower,
                                upper=args.upper,res=args.resolution,
                                width=args.width)
        elif args.balance:
            X = scoreUtils.Chromosome(straw.straw('KR',args.path,cname,cname,'BP',args.resolution),
                                model=model,
                                cname=cname,lower=args.lower,
                                upper=args.upper,res=args.resolution,
                                width=args.width)
        else:
            X = scoreUtils.Chromosome(straw.straw('NONE',args.path,cname,cname,'BP',args.resolution),
                                model=model,
                                cname=cname,lower=args.lower,
                                upper=args.upper,res=args.resolution,
                                width=args.width)
 
        result,R = X.score()
        X.writeBed(args.output,result,R)

