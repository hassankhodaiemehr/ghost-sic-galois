# metric group G = Z^2/Z^2(gamma - I) (row vectors), RW eq. (3)
import numpy as np, json, sys, itertools
from math import gcd
g=np.array(json.loads(sys.argv[1])); N=int(round(np.trace(g)))-2
A=(g-np.eye(2,dtype=int)).astype(int)   # lattice = rows of A (integer combos)
# HNF of rows
def hnf(M):
    M=[list(map(int,r)) for r in M]
    a,b=M[0]; c,d=M[1]
    # make lower-left zero via gcd on first column
    import math
    def ext(x,y):
        if y==0: return (1,0,x) if x>=0 else (-1,0,-x)
        q,r=divmod(x,y); s,t,gg=ext(y,r); return (t,s-q*t,gg)
    s,t,gg=ext(a,c); r1=[s*a+t*c, s*b+t*d]; r2=[(-c//gg)*a+(a//gg)*c, (-c//gg)*b+(a//gg)*d]
    if r2[1]<0: r2=[-r2[0],-r2[1]]
    return r1,r2
r1,r2=hnf(A); assert r1[0]*r2[1]==N or r1[0]*r2[1]==-N, (r1,r2,N)
def red(u):
    u=list(u); k=u[0]//r1[0]; u=[u[0]-k*r1[0],u[1]-k*r1[1]]; u[1]%=abs(r2[1]); return tuple(u)
els=sorted({red((x,y)) for x in range(abs(r1[0])) for y in range(abs(r2[1]))})
assert len(els)==N
idx={e:i for i,e in enumerate(els)}
w=lambda v,u: v[0]*u[1]-v[1]*u[0]
def uG(u): return (np.array(u)@g).tolist()
met=[]; 
for u in els:
    ug=uG(u); sgn=(-1)**((u[0]+u[1]+u[0]*u[1])%2)
    met.append([sgn, w(ug,u)%(2*N)])   # <u> = sgn * e_{2N}(w)
bic=[[w([uG(u)[0]-u[0],uG(u)[1]-u[1]],v)%N for v in els] for u in els]   # <u;v>=e_N(.)
neg=[idx[red((-u[0],-u[1]))] for u in els]
add=[[idx[red((u[0]+v[0],u[1]+v[1]))] for v in els] for u in els]
json.dump(dict(N=N,met=met,bic=bic,neg=neg,add=add),open(sys.argv[2],'w'))
print('N',N,'elements',len(els))
