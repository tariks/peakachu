#!/usr/bin/env python
import pathlib

def main(args):
    import argparse, sys
    import numpy as np
    np.seterr(divide='ignore',invalid='ignore')
    import cooler
    Lib = cooler.Cooler(args.path)
    #resolution = Lib.binsize
    resolution = 10000
    sums=0
    for k in Lib.chromnames[:]:
        x=Lib.matrix(balance=False,sparse=False).fetch(k)
        sums+=np.sum(np.triu(x))

    print('processed {} counts'.format(sums))


