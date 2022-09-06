#!/usr/bin/env python

def main(args):

    import gc, os, joblib
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
    coords = trainUtils.parsebed(args.bedpe, lower=(args.width+1)*res)
    kde, lower, long_start, long_end = trainUtils.learn_distri_kde(coords, res=res)

    if not hic:
        import cooler
        Lib = cooler.Cooler(args.path)
        chromosomes = Lib.chromnames[:]
    else:
        chromosomes = utils.get_hic_chromosomes(args.path, res)

    # train model per chromosome
    collect = {}
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
            neg_coords = trainUtils.negative_generating(
                X, kde, clist, lower, long_start, long_end)
            pos_set = trainUtils.buildmatrix(X, clist, w=args.width)
            neg_set = trainUtils.buildmatrix(X, neg_coords, w=args.width)
            if (not pos_set is None) and (not neg_set is None):
                neg_set = neg_set[:len(pos_set)]
                trainset = pos_set + neg_set
                trainset = np.r_[trainset]
                labels = [1] * len(pos_set) + [0] * len(neg_set)
                labels = np.r_[labels]
                collect[chromname] = [trainset, labels]
            else:
                print(chromname, ' failed to gather fts')
        except:
            print(chromname, ' failed to gather fts')

    for key in chromosomes:
        if key.startswith('chr'):
            chromname = key
        else:
            chromname = 'chr'+key
        
        if chromname != 'chrY':
            continue

        trainset = []
        labels_ = np.r_[[]]
        for ci in collect:
            if (ci != chromname) and (len(collect[ci][1]) > 1):
                trainset.append(collect[ci][0])
                labels_ = np.r_[labels_, collect[ci][1]]
        trainset = np.vstack(trainset)
        
        pn = np.count_nonzero(labels_)
        nn = labels_.size - pn

        print(chromname, 'pos/neg: ', pn, nn)
        model = trainUtils.trainRF(trainset, labels_, nproc=args.nproc)

        joblib.dump(model, os.path.join(args.output, '{0}.pkl'.format(chromname)), compress=('xz', 3))
