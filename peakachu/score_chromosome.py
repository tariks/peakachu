#!/usr/bin/env python

def main(args):

    import argparse, gc, pathlib, os
    import numpy as np
    from sklearn.externals import joblib
    from peakachu import scoreUtils, utils

    np.seterr(divide='ignore',invalid='ignore')
    
    pathlib.Path(args.output).mkdir(parents=True, exist_ok=True) 

    model = joblib.load(args.model)

    hic_info = utils.read_hic_header(args.path) # more robust to check if a file is .hic
    if hic_info is None:
        hic=False
    else:
        hic=True

    if not hic:
        import cooler
        Lib = cooler.Cooler(args.path)
        chromosomes = Lib.chromnames[:]
    else:
        chromosomes = utils.get_hic_chromosomes(args.path, args.resolution)
    
    pre = utils.find_chrom_pre(chromosomes)
    tmp = os.path.split(args.model)[1] # support full path
    ccname = pre + tmp.split('.pk')[0].lstrip('chr') # ccname is consistent with chromosome labels in .hic / .cool
    cikada = 'chr' + ccname.lstrip('chr') # cikada always has prefix "chr"
    
    if not hic:
        X = scoreUtils.Chromosome(Lib.matrix(balance=args.balance, sparse=True).fetch(ccname).tocsr(),
                                model=model,
                                cname=cikada,lower=args.lower,
                                upper=args.upper,res=args.resolution,
                                width=args.width)
    else:
        if args.balance:
            X = scoreUtils.Chromosome(utils.csr_contact_matrix('KR',args.path,ccname,ccname,'BP',args.resolution),
                                model=model,
                                cname=cikada,lower=args.lower,
                                upper=args.upper,res=args.resolution,
                                width=args.width)
        else:
            X = scoreUtils.Chromosome(utils.csr_contact_matrix('NONE',args.path,ccname,ccname,'BP',args.resolution),
                                model=model,
                                cname=cikada,lower=args.lower,
                                upper=args.upper,res=args.resolution,
                                width=args.width)
    result,R = X.score()
    X.writeBed(args.output,result,R)

