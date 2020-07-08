import sys, os

def parse_probs(cell, chrom):

    fil1 = os.path.join('{0}-CTCF'.format(cell), '{0}.bed'.format(chrom))
    fil2 = os.path.join('{0}-H3K27ac'.format(cell), '{0}.bed'.format(chrom))
    D = {}
    with open(fil1, 'r') as source:
        for line in source:
            p = line.rstrip().split()
            loci1 = int(p[1])
            loci2 = int(p[4])
            prob = float(p[6])
            D[(loci1, loci2)] = prob
    
    with open(fil2, 'r') as source:
        for line in source:
            p = line.rstrip().split()
            loci1 = int(p[1])
            loci2 = int(p[4])
            prob = float(p[6])
            if (loci1, loci2) in D:
                D[(loci1, loci2)] = max(D[(loci1, loci2)], prob)
    
    return D

def parse_loops(fil1, fil2):

    loop_set = set()
    with open(fil1, 'r') as source:
        for line in source:
            p = line.rstrip().split()
            loop_set.add(tuple(p[:6]))
    
    with open(fil2, 'r') as source:
        for line in source:
            p = line.rstrip().split()
            loop_set.add(tuple(p[:6]))
    
    return loop_set


chroms = ['chr'+str(i) for i in range(1, 23)] + ['chrX']
infil1 = sys.argv[1] # GM12878
infil2 = sys.argv[2] # another cell
loop_set = parse_loops(infil1, infil2)
pool = []
cell1 = infil1.split('.')[0]
cell2 = infil2.split('.')[0]
for c in chroms:
    print(c)
    probs_1 = parse_probs(cell1, c)
    probs_2 = parse_probs(cell2, c)
    loops = [l for l in loop_set if l[0]==c]
    for l in loops:
        key = (int(l[1]), int(l[4]))
        p1 = probs_1.get(key, 0)
        p2 = probs_2.get(key, 0)
        if p1 or p2:
            pool.append(l+('{0:.4g}'.format(p1), '{0:.4g}'.format(p2)))

with open(sys.argv[3], 'w') as out:
    for line in pool:
        out.write('\t'.join(list(line))+'\n')


