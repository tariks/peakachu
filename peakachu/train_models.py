#!/usr/bin/env python
import pickle
import gc
import pathlib
import random

def main(args):
    import numpy as np
    np.seterr(divide='ignore',invalid='ignore')
    from peakachu import trainUtils
    from scipy.sparse import csr_matrix
    import cooler
    Lib = cooler.Cooler(args.path)
    #resolution = Lib.binsize
    resolution = 10000
    coords = trainUtils.parsebed(args.bedpe,lower=2,res=resolution)
    # train model per chromosome
    positive_class = {}
    negative_class = {}
    for key in Lib.chromnames[:]:
        if key.startswith('chr'):
            chromname=key
        else:
            chromname='chr'+key
        X = Lib.matrix(balance=True,sparse=True).fetch(Key).tocsr() # lower down the memory usage
        R,C = X.nonzero()
        validmask = np.isfinite(X.data) # non-nan
        R,C,data = R[validmask],C[validmask],X.data[validmask]
        X = csr_matrix((data, (R, C)), shape=X.shape)
        del data
        idx = (C-R > 2) & (C-R < 80)
        R,C = R[idx],C[idx]
        idx = np.arange(R.size)
        idx = np.random.choice(idx,100000)
        R,C = R[idx],C[idx]
        clist = coords[chromname]
        try:
            positive_class[chromname] = np.vstack((f for f in trainUtils.buildmatrix(
                                             X,coords[chromname],width=args.width,
                                             chrm=chromname,res=resolution,positive=True)))
            neg_coords = [(r,c) for r,c in zip(R,C)]
            stop = int(chroms[-1].shape[0] * 1.2)
            negative_class[chromname]=np.vstack((f for f in trainUtils.buildmatrix(
                                 X,neg_coords,width=args.width,
                                 chrm=chromname,
                                 res=resolution,positive=False,stop=stop))))

            except:
                print('{} failed with {} coords'.format(Key,len(clist)))
    for key in Lib.chromnames[:]:
        if key.startswith('chr'):
            chromname=key
        else:
            chromname='chr'+key
 
            Xtrain = np.vstack((v for k,v in positive_class.items() if k!=chromname))
            Xfake = np.vstack((v for k,v in negative_class.items() if k!=chromname))
            print(chromname,'pos/neg: ',Xtrain.shape[0],Xfake.shape[0])

            model = trainUtils.trainRF(Xtrain,Xfake)

        with open(args.output+'/'+chromname+'.pkl','wb') as o:
            pickle.dump(model,o)

