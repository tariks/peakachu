#!/usr/env/bin python

def main(args):

    import gc
    import sys
    import numpy as np
    import pandas as pd
    from collections import defaultdict
    from peakachu import peakacluster

    res = args.resolution
    
    clusters, score_pool = peakacluster.parse_peakachu(args.infile, args.threshold, res)
    for c in clusters:
        for p in clusters[c]:
            if p in score_pool[c]:
                s1 = str(p[0]*res)
                e1 = str(p[0]*res+res)
                s2 = str(p[1]*res)
                e2 = str(p[1]*res+res)
                prob = str(score_pool[c][p][0])
                raw_signal = str(score_pool[c][p][1])
                line = [c, s1, e1, c, s2, e2, prob, raw_signal]
                out.write('\t'.join(line)+'\n')
