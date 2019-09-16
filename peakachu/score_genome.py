#!/usr/bin/env python

def main(args):

    import argparse, gc, pathlib
    from sklearn.externals import joblib
    import numpy as np
    from peakachu import scoreUtils, utils

    np.seterr(divide='ignore',invalid='ignore')
    
    pathlib.Path(args.output).mkdir(parents=True, exist_ok=True)

    model = joblib.load(args.model)

    hic_info = utils.read_hic_header(args.path) # more robust to check if a file is .hic
    if hic_info is None:
        hic=False
        import cooler
        Lib = cooler.Cooler(args.path)
        chromosomes = Lib.chromnames[:]
        #nam = args.path.split('.cool')[0]
    else:
        hic=True
        chromosomes = utils.get_hic_chromosomes(args.path, args.resolution)
        #nam = args.path.split('.hic')[0]
    #nam = nam.split('/')[-1]

    for key in chromosomes:
        if key.startswith('chr'):
            cname=key
        else:
            cname='chr'+key
        if not hic:
            X = scoreUtils.Chromosome(Lib.matrix(balance=args.balance, sparse=True).fetch(key).tocsr(),
                                model=model,
                                cname=cname,lower=args.lower,
                                upper=args.upper,res=args.resolution,
                                width=args.width)
        else:
            if args.balance:
                X = scoreUtils.Chromosome(utils.csr_contact_matrix('KR',args.path,key,key,'BP',args.resolution),
                                model=model,
                                cname=cname,lower=args.lower,
                                upper=args.upper,res=args.resolution,
                                width=args.width)
            else:
                X = scoreUtils.Chromosome(utils.csr_contact_matrix('NONE',args.path,key,key,'BP',args.resolution),
                                model=model,
                                cname=cname,lower=args.lower,
                                upper=args.upper,res=args.resolution,
                                width=args.width)
 
        result,R = X.score()
        X.writeBed(args.output,result,R)

