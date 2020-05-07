
def main(args):

    import argparse
    import sys
    import numpy as np
    from scipy import sparse
    from peakachu import utils

    np.seterr(divide='ignore', invalid='ignore')

    # more robust to check if a file is .hic
    check = utils.read_hic_header(args.path)
    if check is None:
        hic = False
    else:
        hic = True

    totals = 0
    if not hic:
        import cooler
        Lib = cooler.Cooler(args.path)

        mindis = args.min_dis // Lib.binsize

        for k in Lib.chromnames[:]:
            intra = np.triu(Lib.matrix(balance=False, sparse=False).fetch(k), k=mindis)
            totals += int(intra.sum())
    else:
        import straw
        hic_info = check
        for k in hic_info['chromsizes']:
            try:
                intra = straw.straw('NONE', args.path, k, k, 'BP', 10000)
                totals += sum(intra[2])  # intra is a list of list x, y, v
            except:
                pass  # handle the inconsistency between .hic header and the matrix

    print(totals)
