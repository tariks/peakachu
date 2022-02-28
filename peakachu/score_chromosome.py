#!/usr/bin/env python

def main(args):

    import pathlib
    import os
    import numpy as np
    from sklearn.externals import joblib
    from peakachu import scoreUtils, utils

    np.seterr(divide='ignore', invalid='ignore')

    pathlib.Path(args.output).mkdir(parents=True, exist_ok=True)

    model = joblib.load(args.model)

    # more robust to check if a file is .hic
    hic_info = utils.read_hic_header(args.path)
    if hic_info is None:
        hic = False
    else:
        hic = True

    if not hic:
        import cooler
        Lib = cooler.Cooler(args.path)

    # ccname is consistent with chromosome labels in .hic / .cool
    ccname = args.chrom
    cikada = 'chr' + ccname.lstrip('chr')  # cikada always has prefix "chr"

    if not hic:
        X = scoreUtils.Chromosome(Lib.matrix(balance=args.balance, sparse=True).fetch(ccname).tocsr(),
                                  model=model,
                                  cname=cikada, lower=args.lower,
                                  upper=args.upper, res=args.resolution,
                                  width=args.width)
    else:
        if args.balance:
            X = scoreUtils.Chromosome(utils.csr_contact_matrix('KR', args.path, ccname, ccname, 'BP', args.resolution),
                                      model=model,
                                      cname=cikada, lower=args.lower,
                                      upper=args.upper, res=args.resolution,
                                      width=args.width)
        else:
            X = scoreUtils.Chromosome(utils.csr_contact_matrix('NONE', args.path, ccname, ccname, 'BP', args.resolution),
                                      model=model,
                                      cname=cikada, lower=args.lower,
                                      upper=args.upper, res=args.resolution,
                                      width=args.width)
    
    result, R = X.score(thre=args.minimum_prob)
    X.writeBed(args.output, result, R)
