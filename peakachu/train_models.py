#!/usr/bin/env python
import pickle
import gc
import pathlib
import random


def main(args):
    import numpy as np
    np.seterr(divide='ignore',invalid='ignore')
    from peakachu import trainUtils
    import cooler
    Lib = cooler.Cooler(args.path)
    #resolution = Lib.binsize
    resolution = 10000
    coords = trainUtils.parsebed(args.bedpe,lower=2,res=resolution)
    chroms = []
    # train model per chromosome
    chromnams=[]
    khromnams=[]
    for k in Lib.chromnames[:]:
        if 'M' not in k and 'Y' not in k:
            khromnams.append(k)
            if k.startswith('c'):
                chromnams.append(k)
            else:
                chromnams.append('chr'+k)
    for k,key in enumerate(khromnams):
        fakechroms=[]
        print('training model for chr{}'.format(key))
        for K,Key in enumerate(khromnams):
            print(key,Key)
            if Key!=key:
                gc.collect()
                X = Lib.matrix(balance=True,sparse=False).fetch(Key)
                R,C = X.nonzero()
                idx = (C-R > 2) & (C-R < 80)
                R,C = R[idx],C[idx]
                idx = np.arange(R.size)
                idx = np.random.choice(idx,100000)
                R,C = R[idx],C[idx]
                np.nan_to_num(X,copy=False)
                #clist = coords['chr'+Key]
                clist = coords[chromnams[K]]
                try:
                    chroms.append(np.vstack((f for f in trainUtils.buildmatrix(
                                             X,coords[chromnams[K]],width=args.width,
                                             chrm=Key,res=resolution,positive=True))))
                    rsize = R.size
                    neg_coords = [(R[j],C[j]) for j in range(rsize)]
                    stop = int(chroms[-1].shape[0] * 1.2)
                    fakechroms.append(np.vstack((f for f in trainUtils.buildmatrix(
                                 X,neg_coords,width=args.width,
                                 chrm=Key,
                                 res=resolution,positive=False,stop=stop))))

                except:
                    print('{} failed with {} coords'.format(Key,len(clist)))
        Xtrain = np.vstack(chroms)
        Xfake = np.vstack(fakechroms)
        del chroms[:]
        del fakechroms[:]
        np.random.shuffle(Xtrain)
        np.random.shuffle(Xfake)
        print('pos/neg: ',Xtrain.shape[0],Xfake.shape[0])

        model = trainUtils.trainRF(Xtrain,Xfake)

        with open(args.output+'/'+chromnams[k]+'.pkl','wb') as o:
            pickle.dump(model,o)

