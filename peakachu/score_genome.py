#!/usr/bin/env python

def main(args):

    import joblib, os
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
        import cooler
        Lib = cooler.Cooler(args.path)
        chromosomes = Lib.chromnames[:]
        #nam = args.path.split('.cool')[0]
    else:
        hic = True
        chromosomes = utils.get_hic_chromosomes(args.path, args.resolution)
        #nam = args.path.split('.hic')[0]
    #nam = nam.split('/')[-1]

    queue = []
    for key in chromosomes:

        chromlabel = key.lstrip('chr')
        if (not args.chroms) or (chromlabel.isdigit() and '#' in args.chroms) or (chromlabel in args.chroms):
            queue.append(key)
    
    for key in queue:

        if key.startswith('chr'):
            cname = key
        else:
            cname = 'chr'+key
            
        if not hic:
            if args.balance:
                M = utils.tocsr(Lib.matrix(balance=args.balance, sparse=True).fetch(key))
                raw_M = utils.tocsr(Lib.matrix(balance=False, sparse=True).fetch(key))
                weights = Lib.bins().fetch(key)['weight'].values
                X = scoreUtils.Chromosome(M, model=model, raw_M=raw_M, weights=weights,
                                          cname=cname, lower=args.lower,
                                          upper=args.upper, res=args.resolution,
                                          width=width)
            else:
                M = utils.tocsr(Lib.matrix(balance=False, sparse=True).fetch(key))
                X = scoreUtils.Chromosome(M, model=model, raw_M=M, weights=None,
                                          cname=cname, lower=args.lower,
                                          upper=args.upper, res=args.resolution,
                                          width=width)
        else:
            if args.balance:
                M = utils.csr_contact_matrix('KR', args.path, key, key, 'BP', args.resolution)
                raw_M = utils.csr_contact_matrix('NONE', args.path, key, key, 'BP', args.resolution)
                X = scoreUtils.Chromosome(M, model=model, raw_M=raw_M, weights=None,
                                          cname=cname, lower=args.lower,
                                          upper=args.upper, res=args.resolution,
                                          width=width)
            else:
                M = utils.csr_contact_matrix('NONE', args.path, key, key, 'BP', args.resolution)
                X = scoreUtils.Chromosome(M, model=model, raw_M=M, weights=None,
                                          cname=cname, lower=args.lower,
                                          upper=args.upper, res=args.resolution,
                                          width=width)

        result, R = X.score(thre=args.minimum_prob)
        X.writeBed(args.output, result, R)
