#!/usr/bin/env python

def main(args):
    
    import os, joblib
    import numpy as np
    from peakachu import scoreUtils, utils

    np.seterr(divide='ignore', invalid='ignore')

    if os.path.exists(args.output):
        os.remove(args.output)

    model = joblib.load(args.model)

    # deduce the width parameter used during the training
    width = int((np.sqrt(model.feature_importances_.size) - 1) / 2)

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
        if args.balance:
            M = utils.tocsr(Lib.matrix(balance=args.balance, sparse=True).fetch(ccname))
            raw_M = utils.tocsr(Lib.matrix(balance=False, sparse=True).fetch(ccname))
            weights = Lib.bins().fetch(ccname)['weight'].values
            X = scoreUtils.Chromosome(M, model=model, raw_M=raw_M, weights=weights,
                                      cname=cikada, lower=args.lower,
                                      upper=args.upper, res=args.resolution,
                                      width=width)
        else:
            M = utils.tocsr(Lib.matrix(balance=False, sparse=True).fetch(ccname))
            X = scoreUtils.Chromosome(M, model=model, raw_M=M, weights=None,
                                      cname=cikada, lower=args.lower,
                                      upper=args.upper, res=args.resolution,
                                      width=width)
    else:
        if args.balance:
            M = utils.csr_contact_matrix('KR', args.path, ccname, ccname, 'BP', args.resolution)
            raw_M = utils.csr_contact_matrix('NONE', args.path, ccname, ccname, 'BP', args.resolution)
            X = scoreUtils.Chromosome(M, model=model, raw_M=raw_M, weights=None,
                                      cname=cikada, lower=args.lower,
                                      upper=args.upper, res=args.resolution,
                                      width=width)
        else:
            M = utils.csr_contact_matrix('NONE', args.path, ccname, ccname, 'BP', args.resolution)
            X = scoreUtils.Chromosome(M, model=model, raw_M=M, weights=None,
                                      cname=cikada, lower=args.lower,
                                      upper=args.upper, res=args.resolution,
                                      width=width)
    
    result, R = X.score(thre=args.minimum_prob)
    X.writeBed(args.output, result, R)
