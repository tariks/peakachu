import sys
import numpy as np

infil = sys.argv[1] # original output from Fit-Hi-C
top_n = int(sys.argv[2]) # an integer
outfil = sys.argv[3] # output file name

qvalues = np.loadtxt(infil, usecols=[-1])
minq = qvalues[qvalues>0].min()

pool = []
with open(infil, 'r') as source:
    for line in source:
        parse = line.rstrip().split()
        count = int(parse[6])
        qvalue = float(parse[7])
        if qvalue <= 0:
            qvalue = minq
        record = (-np.log(qvalue), count) + tuple(parse[:6])
        pool.append(record)

pool.sort(reverse=True)

selected = pool[:top_n]
with open(outfil, 'w') as out:
    for line in selected:
        out.write('\t'.join(list(line[2:])+[str(line[0])])+'\n')