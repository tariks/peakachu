import sys
from collections import defaultdict

bedfile = sys.argv[1]
ctcfsites = 'ctcf.orient'

convergent = 0
convergent2 = 0
singles = 0
F=defaultdict(list)
R=defaultdict(list)
for line in open(ctcfsites):
    s = line.split()
    c,a,orient = s[0],int(s[1]),s[2]
    if orient == 'Forward':
        F[c].append(a)
    if orient == 'Reverse':
        R[c].append(a)
#for k in R.keys():
#    R[k] = list(set(R[k]))
#    F[k] = list(set(F[k]))
total=0
doubles=0
peakloci=defaultdict(list)
for line in open(bedfile):
    total+=1
    s = line.split()
    c,a,b = s[0],int(s[1]),int(s[4])
    a-=a%10000
    b-=b%10000
    peakloci[c]+=[a,b]
    aF,aR,bF,bR = 0,0,0,0
    aF+=F[c].count(a)
    aR+=R[c].count(a)
    bF+=F[c].count(b)
    bR+=R[c].count(b)
    if aF+bF+aR+bR > 0:
        singles+=1
    if aF>=1 and bR>=1:
        convergent2+=1
    if aF+aR==1 and bF+bR==1:
        doubles+=1
        if aF > 0 and bR > 0:
            convergent+=1
for k,v in peakloci.items():
    peakloci[k] = list(set(v))
ptotal,pnum=0,0
for k,v in peakloci.items():
    ptotal+=len(v)
    for i in v:
        if i in F[k] or i in R[k]:
            pnum+=1
print(total,'total sites')
print(doubles,'doubled sites')
if doubles > 0:
    print(convergent/doubles,'convergent ratio')
print(convergent2,' all convergent peaks')
print(pnum/ptotal,' CTCF loci')
print(singles/total,' min 1 ctcf')

