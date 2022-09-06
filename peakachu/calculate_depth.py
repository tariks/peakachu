
def main(args):

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
        genome_size = Lib.chromsizes.sum()

        mindis = args.min_dis // Lib.binsize

        for k in Lib.chromnames[:]:
            print(k)
            intra = np.triu(Lib.matrix(balance=False, sparse=False).fetch(k), k=mindis)
            totals += int(intra.sum())
        
    else:
        import straw
        hic_info = check
        genome_size = sum([hic_info['chromsizes'][c] for c in hic_info['chromsizes']])

        for k in hic_info['chromsizes']:
            print(k)
            try:
                lowres = max(hic_info['Base pair-delimited resolutions'])
                intra = straw.straw('NONE', args.path, k, k, 'BP', lowres)
                totals += sum(intra[2])  # intra is a list of list x, y, v
            except:
                print('chrom', k, 'failed')
                pass  # handle the inconsistency between .hic header and the matrix

    print('num of intra reads in your data:', totals)
    matched_read_num = 3031042417 / genome_size * totals
    print('num of intra reads in a human with matched sequencing coverage:', int(matched_read_num))
    print('suggested model:', match_pretrained_models(matched_read_num))

def match_pretrained_models(v, platform='Hi-C'):

    import numpy as np

    if platform in ['Hi-C', 'Micro-C']:
        arr = [
            5000000, 10000000, 30000000, 50000000, 100000000,
            150000000, 200000000, 250000000, 300000000, 350000000,
            400000000, 450000000, 500000000, 550000000, 600000000,
            650000000, 700000000, 750000000, 800000000, 850000000,
            900000000, 1000000000, 1200000000, 1400000000, 1600000000,
            1800000000, 2000000000
        ]

    diff = np.abs(v - np.r_[arr])
    idx = np.argmin(diff)
    if arr[idx] >= 1000000000:
        label = '{0:.2g} billion'.format(arr[idx]/1000000000)
    else:
        label = '{0} million'.format(arr[idx]//1000000)

    return label