import numpy as np, itertools, json, subprocess
S=np.array([[0,-1],[1,0]]); T=lambda b: np.array([[1,b],[0,1]])
seen=set(); out=[]
def canon(c):
    rots=[c[i:]+c[:i] for i in range(len(c))]; return min(rots)
for L in range(1,5):
    for c in itertools.product(range(2,40),repeat=L):
        if all(x==2 for x in c): continue
        cc=canon(list(c))
        if tuple(cc) in seen: continue
        g=np.eye(2,dtype=int)
        for b in cc: g=g@T(b)@S
        N=int(np.trace(g))-2
        if N<5 or N>30: continue
        seen.add(tuple(cc)); out.append((N,cc,g.tolist()))
out.sort()
print(len(out))
json.dump(out,open('gammas.json','w'))
