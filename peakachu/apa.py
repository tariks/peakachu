#!/usr/bin/env python
import cooler, argparse, sys, logging, logging.handlers
import numpy as np
import pandas as pd
import random
import sys
#currentVersion = peakachu.__version__
currentVersion = '0.2'
np.seterr(divide='ignore',invalid='ignore')

p=sys.argv[1]
loops = sys.argv[2]

import cooler
Lib = cooler.Cooler(p)
resolution = 10000
x=pd.read_table(loops,index_col=0,usecols=[0,1,4])
chroms = list(set(x.index))
apa=np.zeros((21,21))
for k in chroms:
    print(k)
    X=x.loc[k].values
    X=X.astype(int)
    X//=10000
    r,c = X[:,0],X[:,1]
    X = Lib.matrix(balance=True,sparse=False).fetch(k.split('chr')[1])
    np.nan_to_num(X,copy=False)
    for i in range(r.size):
        try:
            apa+=X[r[i]-10:r[i]+11,c[i]-10:c[i]+11]
        except:
            pass
np.savetxt(sys.argv[3],apa)

