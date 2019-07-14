from collections import Counter
import numpy as np
from scipy.signal import find_peaks, peak_widths
from sklearn.cluster import dbscan
from scipy.spatial.distance import euclidean

def find_anchors(pos, min_count=3, min_dis=20000, wlen=800000, res=10000):

    min_dis = max(min_dis//res, 1)
    wlen = min(wlen//res, 20)

    count = Counter(pos)
    refidx = range(min(count), max(count)+1)
    signal = np.r_[[count[i] for i in refidx]]
    summits = find_peaks(signal, height=min_count, distance=min_dis)[0]
    sorted_summits = [(signal[i],i) for i in summits]
    sorted_summits.sort(reverse=True) # sort by peak count
    
    peaks = set()
    records = {}
    for _, i in sorted_summits:
        tmp = peak_widths(signal, [i], rel_height=1, wlen=wlen)[2:4]
        li, ri = int(np.round(tmp[0][0])), int(np.round(tmp[1][0]))
        lb = refidx[li]
        rb = refidx[ri]
        if not len(peaks):
            peaks.add((refidx[i], lb, rb))
            for b in range(lb, rb+1):
                records[b] = (refidx[i], lb, rb)
        else:
            for b in range(lb, rb+1):
                if b in records:
                    # merge anchors
                    m_lb = min(lb, records[b][1])
                    m_rb = max(rb, records[b][2])
                    summit = records[b][0] # always the highest summit
                    peaks.remove(records[b])
                    break
            else: # loop terminates normally
                m_lb, m_rb, summit = lb, rb, refidx[i]
            peaks.add((summit, m_lb, m_rb))
            for b in range(m_lb, m_rb+1):
                records[b] = (summit, m_lb, m_rb)
    
    return peaks

def _cluster_core(sort_list, r, visited, final_list):

    pos = np.r_[[i[1] for i in sort_list]]
    if len(pos) >= 2:
        _, labels = dbscan(pos, eps=r, min_samples=2)
        pool = set()
        for i, p in enumerate(sort_list):
            if p[1] in pool:
                continue
            c = labels[i]
            if c==-1:
                continue
            sub = pos[labels==c]
            cen = p[1]
            rad = r
            Local = [p[1]]
            ini = -1
            while len(sub):
                out = []
                for q in sub:
                    if tuple(q) in pool:
                        continue
                    tmp = euclidean(q, cen)
                    if tmp<=rad:
                        Local.append(tuple(q))
                    else:
                        out.append(tuple(q))
                if len(out)==ini:
                    break
                ini = len(out)
                tmp = np.r_[Local]
                # assign centroid to a certain pixel
                cen = tuple(tmp.mean(axis=0).round().astype(int))
                rad = np.int(np.round(max([euclidean(cen,q) for q in Local]))) + r
                sub = np.r_[out]
            for q in Local:
                pool.add(q)
            final_list.append((p[1], cen, rad))
        
        visited.update(pool)

def local_clustering(Donuts, res, min_count=3, r=20000, sumq=1):

    final_list = []
    x = np.r_[[i[0] for i in Donuts]]
    y = np.r_[[i[1] for i in Donuts]]
    if x.size == 0:
        return final_list

    x_anchors = find_anchors(x, min_count=min_count, min_dis=r, res=res)
    y_anchors = find_anchors(y, min_count=min_count, min_dis=r, res=res)
    r = max(r//res, 1)
    visited = set()
    lookup = set(zip(x, y))
    for x_a in x_anchors:
        for y_a in y_anchors:
            sort_list = []
            for i in range(x_a[1], x_a[2]+1):
                for j in range(y_a[1], y_a[2]+1):
                    if (i, j) in lookup:
                        sort_list.append((Donuts[(i,j)], (i,j)))
            sort_list.sort(reverse=True)
            _cluster_core(sort_list, r, visited, final_list)
    
    sort_list = [] # out of anchor
    for i, j in zip(x, y):
        if (i,j) in visited:
            continue
        sort_list.append((Donuts[(i,j)], (i,j)))
    sort_list.sort(reverse=True)
    _cluster_core(sort_list, r, visited, final_list)

    x_summits = set([i[0] for i in x_anchors])
    y_summits = set([i[0] for i in y_anchors])
    for i, j in zip(x, y):
        if (i,j) in visited:
            continue
        
        if (i in x_summits) or (j in y_summits):
            final_list.append(((i,j), (i,j), 0))
    
    return final_list