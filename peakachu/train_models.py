#!/usr/bin/env python

def main(args):

    from sklearn.externals import joblib
    import gc
    import pathlib
    import straw
    import numpy as np
    from peakachu import trainUtils, utils

    np.seterr(divide='ignore', invalid='ignore')

    pathlib.Path(args.output).mkdir(parents=True, exist_ok=True)

    # more robust to check if a file is .hic
    hic_info = utils.read_hic_header(args.path)

    if hic_info is None:
        hic = False
    else:
        hic = True

    res = args.resolution
    coords = trainUtils.parsebed(args.bedpe, lower=6*res)
    kde, lower, long_start, long_end = trainUtils.learn_distri_kde(coords, res=res)

    if not hic:
        import cooler
        Lib = cooler.Cooler(args.path)
        chromosomes = Lib.chromnames[:]
    else:
        chromosomes = utils.get_hic_chromosomes(args.path, res)

    # train model per chromosome
    positive_class = {}
    negative_class = {}
    for key in chromosomes:
        if key.startswith('chr'):
            chromname = key
        else:
            chromname = 'chr'+key
        print('collecting from {}'.format(key))
        if not hic:
            X = Lib.matrix(balance=args.balance,
                           sparse=True).fetch(key).tocsr()
        else:
            if args.balance:
                X = utils.csr_contact_matrix(
                    'KR', args.path, key, key, 'BP', res)
            else:
                X = utils.csr_contact_matrix(
                    'NONE', args.path, key, key, 'BP', res)

        # deal with the situation when resolutions of the matrix
        # and the training set are different
        clist = []
        for s1, e1, s2, e2 in coords[chromname]:
            bins1 = range(s1//res, (e1+res-1)//res)
            bins2 = range(s2//res, (e2+res-1)//res)
            maxv = 0
            binpair = None
            for b1 in bins1:
                for b2 in bins2:
                    if X[b1, b2] > maxv:
                        maxv = X[b1, b2]
                        binpair = (b1, b2)
            if maxv > 0:
                clist.append(binpair)

        try:
            positive_class[chromname] = np.vstack((f for f in trainUtils.buildmatrix(
                X, clist, width=args.width)))
            neg_coords = trainUtils.negative_generating(
                X, kde, clist, lower, long_start, long_end)
            stop = len(clist)
            negative_class[chromname] = np.vstack((f for f in trainUtils.buildmatrix(
                X, neg_coords, width=args.width,
                positive=False, stop=stop)))
        except:
            print(chromname, ' failed to gather fts')

    for key in chromosomes:
        if key.startswith('chr'):
            chromname = key
        else:
            chromname = 'chr'+key

        Xtrain = np.vstack(
            (v for k, v in positive_class.items() if k != chromname))
        Xfake = np.vstack(
            (v for k, v in negative_class.items() if k != chromname))
        print(chromname, 'pos/neg: ', Xtrain.shape[0], Xfake.shape[0])
        model = trainUtils.trainRF(Xtrain, Xfake)

        joblib.dump(model, args.output+'/'+chromname +
                    '.pkl', compress=('xz', 3))
