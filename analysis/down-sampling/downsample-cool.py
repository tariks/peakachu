import cooler, random, sys, os

cool_uri = sys.argv[1]
lib = cooler.Cooler(cool_uri)
rate = float(sys.argv[2]) # specify the down-sampling rate here, 0.03, 0.1, ...
outdir = sys.argv[3] # the output folder name

if not os.path.exists(outdir):
    os.mkdir(outdir)

chroms = ['chr'+str(i) for i in range(1, 23)] + ['chrX']

for c in chroms:
    M = lib.matrix(balance=False, sparse=True).fetch(c)
    x, y = M.nonzero()
    data = M.data
    outfil = os.path.join(outdir, '{0}_{0}.txt'.format(c.lstrip('chr')))
    with open(outfil, 'w') as out:
        for i, j, v in zip(x, y, data):
            if i > j:
                continue
            count = 0
            for t in range(int(v)):
                rnd = random.random()
                if rnd < rate:
                    count += 1
            if count > 0:
                line = [str(i), str(j), str(count)]
                out.write('\t'.join(line)+'\n')
        