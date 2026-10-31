# Zauner-symmetric U_P-ghost search for odd d: x in an eigenspace of the Clifford unitary of F=[[0,-1],[1,-1]].
import numpy as np, sys
from scipy.optimize import least_squares
d=int(sys.argv[1]); seed=int(sys.argv[2]) if len(sys.argv)>2 else 0
w=np.exp(2j*np.pi/d); tau=-np.exp(1j*np.pi/d)
X=np.roll(np.eye(d),1,axis=0); Z=np.diag(w**np.arange(d)); UP=np.eye(d)[[(-j)%d for j in range(d)]]
def Dp(a,b): return (tau**((a*b)%d))*np.linalg.matrix_power(X,a%d)@np.linalg.matrix_power(Z,b%d)
al,be,ga,de=0,-1,1,-1   # F
bi=pow(be%d,-1,d)
U=np.array([[tau**((bi*(al*s*s-2*r*s+de*r*r))%d) for s in range(d)] for r in range(d)])/np.sqrt(d)
# check covariance U D_p U^dag = D_{Fp} up to phase
for (a,b) in [(1,0),(0,1),(1,1)]:
    Fp=((al*a+be*b)%d,(ga*a+de*b)%d)
    M=U@Dp(a,b)@U.conj().T; R=Dp(*Fp)
    ph=np.vdot(R.flatten(),M.flatten())/d
    assert abs(abs(ph)-1)<1e-9 and np.linalg.norm(M-ph*R)<1e-8, ("cov",a,b,ph)
ev,V=np.linalg.eig(U)
ps=[(a,b) for a in range(d) for b in range(d) if (a,b)!=(0,0)]
neg={p:((-p[0])%d,(-p[1])%d) for p in ps}
half=[p for p in ps if ps.index(p)<ps.index(neg[p])]
B={p:UP@Dp(*p) for p in ps}
# group eigenvalues
cls=[]
for k in range(d):
    for c in cls:
        if abs(ev[c[0]]-ev[k])<1e-6: c.append(k); break
    else: cls.append([k])
rng=np.random.default_rng(seed)
for trial in range(4000):
    c=cls[trial%len(cls)]; Q,_=np.linalg.qr(V[:,c]); m=len(c)
    def res(v):
        x=Q@(v[:m]+1j*v[m:]); n=np.vdot(x,UP@x)
        r=[]
        for p in half:
            q=(d+1)*np.vdot(x,B[p]@x)*np.vdot(x,B[neg[p]]@x)/n**2-1; r+=[q.real,q.imag]
        r+=[np.vdot(x,x).real-1]
        return np.array(r)
    v0=rng.normal(size=2*m)
    try: r=least_squares(res,v0,xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=3000)
    except Exception: continue
    if np.max(np.abs(r.fun))<1e-10:
        x=Q@(r.x[:m]+1j*r.x[m:]); n=np.vdot(x,UP@x); Pi=np.outer(x,(UP@x).conj())/n
        herm=np.linalg.norm(Pi-Pi.conj().T)
        print('trial',trial,'eigsp dim',m,'res',np.max(np.abs(r.fun)),'non-herm',herm,file=sys.stderr,flush=True)
        if herm>1e-3:
            with open(f'ghost{d}.txt','w') as fo:
                fo.write(' '.join(repr(float(z.real)) for z in x)+'\n'+' '.join(repr(float(z.imag)) for z in x)+'\n')
            break
