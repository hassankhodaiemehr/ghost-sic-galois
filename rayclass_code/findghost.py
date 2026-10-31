import numpy as np, sys
from scipy.optimize import least_squares
d=int(sys.argv[1]); seed=int(sys.argv[2]) if len(sys.argv)>2 else 0
w=np.exp(2j*np.pi/d); tau=-np.exp(1j*np.pi/d)
X=np.roll(np.eye(d),1,axis=0); Z=np.diag(w**np.arange(d)); UP=np.eye(d)[[(-j)%d for j in range(d)]]
ps=[(a,b) for a in range(d) for b in range(d) if (a,b)!=(0,0)]
D={p:(tau**(p[0]*p[1]))*np.linalg.matrix_power(X,p[0])@np.linalg.matrix_power(Z,p[1]) for p in ps}
neg={p:((-p[0])%d,(-p[1])%d) for p in ps}
half=[p for p in ps if ps.index(p)<ps.index(neg[p]) or p==neg[p]]
B={p:UP@D[p] for p in ps}
def res(v):
    x=v[:d]+1j*v[d:]; n=np.vdot(x,UP@x)
    a={p:np.vdot(x,B[p]@x)/n for p in ps}
    r=[]
    for p in half:
        q=(d+1)*a[p]*a[neg[p]]-1; r+=[q.real,q.imag]
    r+=[ (np.vdot(x,x).real-1) ]
    return np.array(r)
rng=np.random.default_rng(seed)
for trial in range(2000):
    v0=rng.normal(size=2*d)
    try:
        r=least_squares(res,v0,xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
    except Exception as ex: continue
    if np.max(np.abs(r.fun))<1e-11:
        x=r.x[:d]+1j*r.x[d:]
        # ghost (not SIC): check non-Hermitian
        n=np.vdot(x,UP@x); Pi=np.outer(x,(UP@x).conj())/n  # Pi = x x^dag U_P / n  (as x y^T with y = conj(U_P x))
        herm=np.linalg.norm(Pi-Pi.conj().T)
        print('trial',trial,'res',np.max(np.abs(r.fun)),'non-herm',herm,file=sys.stderr)
        if herm>1e-3:
            with open(f'ghost{d}.txt','w') as fo:
                fo.write(' '.join(repr(float(z.real)) for z in x)+'\n'+' '.join(repr(float(z.imag)) for z in x)+'\n')
            break
