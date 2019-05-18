#!/usr/bin/env python
from sklearn.externals import joblib
import gc
import pathlib

def main(args):
    import numpy as np
    np.seterr(divide='ignore',invalid='ignore')
    from peakachu import trainUtils
    from scipy.sparse import csr_matrix
    import cooler
    import random
    pathlib.Path(args.output).mkdir(parents=True, exist_ok=True)
    Lib = cooler.Cooler(args.path)
    #resolution = Lib.binsize
    resolution = args.resolution
    coords = trainUtils.parsebed(args.bedpe,lower=1,res=resolution)
    negcoords = trainUtils.parsebed(args.control,lower=1,res=resolution)
    # train model per chromosome
    positive_class = {}
    negative_class = {}
    for key in Lib.chromnames[:]:
        if key.startswith('chr'):
            chromname=key
        else:
            chromname='chr'+key
        print('collecting from {}'.format(key))
        X = Lib.matrix(balance=args.balance,sparse=True).fetch(key).tocsr() # lower down the memory usage
        R,C = X.nonzero()
        validmask = np.isfinite(X.data) # non-nan
        R,C,data = R[validmask],C[validmask],X.data[validmask]
        X = csr_matrix((data, (R, C)), shape=X.shape)
        del data
        clist = coords[chromname]
        nlist = negcoords[chromname]
        try:
            positive_class[chromname] = np.vstack((f for f in trainUtils.buildmatrix(
                                             X,clist,width=args.width)))
            stop=positive_class[chromname].shape[0]
            negative_class[chromname] = np.vstack((f for f in trainUtils.buildmatrix(
                                             X,nlist,stop=stop,positive=False,width=args.width)))
        except:
            print(chromname, ' failed to gather fts')

    for key in Lib.chromnames[:]:
        if key.startswith('chr'):
            chromname=key
        else:
            chromname='chr'+key
 
        Xtrain = np.vstack((v for k,v in positive_class.items() if k!=chromname))
        Xfake = np.vstack((v for k,v in negative_class.items() if k!=chromname))
        print(chromname,'pos/neg: ',Xtrain.shape[0],Xfake.shape[0])
        model = trainUtils.trainRF(Xtrain,Xfake)

        joblib.dump(model, args.output+'/'+chromname+'.pkl', compress=('xz',3))

