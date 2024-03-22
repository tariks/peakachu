
# Read information from the hic header

import struct, io, os
import numpy as np
from sklearn.isotonic import IsotonicRegression
from scipy.sparse import csr_matrix
from numba import njit

def tocsr(X):

    row, col, data = X.row, X.col, X.data
    M = csr_matrix((data, (row, col)), shape=X.shape, dtype=float)

    return M

def csr_contact_matrix(norm, hicfile, chr1loc, chr2loc, unit,
                       binsize, is_synapse=False):
    '''
    Extract the contact matrix from .hic in CSR sparse format
    '''
    import straw

    tri_list = straw.straw(norm, hicfile, chr1loc,
                           chr2loc, unit, binsize, is_synapse)
    row = [r//binsize for r in tri_list[0]]
    col = [c//binsize for c in tri_list[1]]
    value = tri_list[2]
    N = max(col) + 1

    # re-scale KR matrix to ICE-matrix range
    M = csr_matrix((value, (row, col)), shape=(N, N), dtype=float)
    margs = np.array(M.sum(axis=0)).ravel() + \
        np.array(M.sum(axis=1)).ravel() - M.diagonal(0)
    margs[np.isnan(margs)] = 0
    scale = margs[margs != 0].mean()
    row, col = M.nonzero()
    value = M.data / scale
    M = csr_matrix((value, (row, col)), shape=(N, N), dtype=float)

    return M

def get_hic_chromosomes(hicfile, res):

    import straw
    
    hic_info = read_hic_header(hicfile)
    chromosomes = []
    # handle with inconsistency between .hic header and matrix data
    for c, Len in hic_info['chromsizes'].items():
        try:
            loc = '{0}:{1}:{2}'.format(c, 0, min(Len, 100000))
            _ = straw.straw('NONE', hicfile, loc, loc, 'BP', res)
            chromosomes.append(c)
        except:
            pass

    return chromosomes


def find_chrom_pre(chromlabels):

    ini = chromlabels[0]
    if ini.startswith('chr'):
        return 'chr'

    else:
        return ''


def readcstr(f):
    buf = ""
    while True:
        b = f.read(1)
        b = b.decode('utf-8', 'backslashreplace')
        if b is None or b == '\0':
            return str(buf)
        else:
            buf = buf + b


def read_hic_header(hicfile):

    if not os.path.exists(hicfile):
        return None  # probably a cool URI

    req = open(hicfile, 'rb')
    magic_string = struct.unpack('<3s', req.read(3))[0]
    req.read(1)
    if (magic_string != b"HIC"):
        return None  # this is not a valid .hic file

    info = {}
    version = struct.unpack('<i', req.read(4))[0]
    info['version'] = str(version)

    masterindex = struct.unpack('<q', req.read(8))[0]
    info['Master index'] = str(masterindex)

    genome = ""
    c = req.read(1).decode("utf-8")
    while (c != '\0'):
        genome += c
        c = req.read(1).decode("utf-8")
    info['Genome ID'] = str(genome)

    nattributes = struct.unpack('<i', req.read(4))[0]
    attrs = {}
    for i in range(nattributes):
        key = readcstr(req)
        value = readcstr(req)
        attrs[key] = value
    info['Attributes'] = attrs

    nChrs = struct.unpack('<i', req.read(4))[0]
    chromsizes = {}
    for i in range(nChrs):
        name = readcstr(req)
        length = struct.unpack('<i', req.read(4))[0]
        if name != 'ALL':
            chromsizes[name] = length

    info['chromsizes'] = chromsizes

    info['Base pair-delimited resolutions'] = []
    nBpRes = struct.unpack('<i', req.read(4))[0]
    for i in range(nBpRes):
        res = struct.unpack('<i', req.read(4))[0]
        info['Base pair-delimited resolutions'].append(res)

    info['Fragment-delimited resolutions'] = []
    nFrag = struct.unpack('<i', req.read(4))[0]
    for i in range(nFrag):
        res = struct.unpack('<i', req.read(4))[0]
        info['Fragment-delimited resolutions'].append(res)

    return info

def calculate_expected(M, maxdis, raw=False):

    n = M.shape[0]
    R, C = M.nonzero()
    valid_pixels = np.isfinite(M.data)
    # extract valid columns
    if raw:
        R, C, data = R[valid_pixels], C[valid_pixels], M.data[valid_pixels]
        M = csr_matrix((data, (R, C)), shape=M.shape, dtype=float)
        marg = np.array(M.sum(axis=0)).ravel()
        valid_cols = marg > 0
    else:
        R, C = set(R[valid_pixels]), set(C[valid_pixels])
        valid_cols = np.zeros(n, dtype=bool)
        for i in R:
            valid_cols[i] = True
        for i in C:
            valid_cols[i] = True
    
    # calculate the expected value for each genomic distance
    exp_arr = np.zeros(maxdis+1)
    for i in range(maxdis+1):
        if i == 0:
            valid = valid_cols
        else:
            valid = valid_cols[:-i] * valid_cols[i:]
        
        diag = M.diagonal(i)
        diag = diag[valid]
        if diag.size > 10:
            exp = diag.mean()
            exp_arr[i] = exp
    
    # make exp_arr stringently non-increasing
    IR = IsotonicRegression(increasing=False, out_of_bounds='clip')
    _d = np.where(exp_arr > 0)[0]
    IR.fit(_d, exp_arr[_d])
    exp_arr = IR.predict(list(range(maxdis+1)))

    return exp_arr

@njit
def distance_normaize_core(sub, exp_bychrom, x, y, w):

    # calculate x and y indices
    x_arr = np.arange(x-w, x+w+1).reshape((2*w+1, 1))
    y_arr = np.arange(y-w, y+w+1)

    D = y_arr - x_arr
    D = np.abs(D)
    min_dis = D.min()
    max_dis = D.max()
    if max_dis >= exp_bychrom.size:
        return sub
    else:
        exp_sub = np.zeros(sub.shape)
        for d in range(min_dis, max_dis+1):
            xi, yi = np.where(D==d)
            for i, j in zip(xi, yi):
                exp_sub[i, j] = exp_bychrom[d]
            
        normed = sub / exp_sub

        return normed

@njit
def image_normalize(arr_2d):

    arr_2d = (arr_2d - arr_2d.min()) / (arr_2d.max() - arr_2d.min()) # value range: [0,1]

    return arr_2d

@njit
def distance_normalize(arr_pool, exp_bychrom, xi, yi, w):

    clist = []
    fea = []
    for i in range(xi.size):
        x = xi[i]
        y = yi[i]
        window = arr_pool[i]

        bad_x, bad_y = np.where(np.isnan(window))
        for i_, j_ in zip(bad_x, bad_y):
            window[i_, j_] = 0
             
        if np.count_nonzero(window) < window.size*0.1:
            continue
        
        ll_mean = window[:w, :w].mean()
        if ll_mean > 0:
            center = window[w, w]
            p2LL = center / ll_mean
            if p2LL > 0.1:
                window = distance_normaize_core(window, exp_bychrom, x, y, w)
                fea.append(window)
                clist.append((x, y))
    
    return fea, clist