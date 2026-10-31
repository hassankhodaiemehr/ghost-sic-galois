import numpy as np, sys
from scipy.optimize import minimize
d=int(sys.argv[1]); seed=int(sys.argv[2]) if len(sys.argv)>2 else 0
w=np.exp(2j*np.pi/d); tau=-np.exp(1j*np.pi/d)
X=np.roll(np.eye(d),1,axis=0); Z=np.diag(w**np.arange(d))
Ds=[(tau**(a*b))*np.linalg.matrix_power(X,a)@np.linalg.matrix_power(Z,b) for a in range(d) for b in range(d)]
A=np.array(Ds)
def f(v):
    psi=v[:d]+1j*v[d:]; n=np.vdot(psi,psi).real
    ov=np.einsum('i,pij,j->p',psi.conj(),A,psi)
    val=np.sum(np.abs(ov)**4)/n**4
    # gradient
    g=np.zeros(d,complex)
    for p in range(d*d):
        g+=2*np.abs(ov[p])**2*(np.conj(ov[p])*(A[p]@psi)+ov[p]*(A[p].conj().T@psi))
    g=g/n**4 - 4*val*psi/n
    return val, np.concatenate([2*g.real,2*g.imag])
target=2*d/(d+1)
rng=np.random.default_rng(seed)
for trial in range(400):
    v0=rng.normal(size=2*d)
    r=minimize(f,v0,jac=True,method='L-BFGS-B',options={'maxiter':5000,'gtol':1e-14,'ftol':1e-16})
    if abs(r.fun-target)<1e-11:
        psi=r.x[:d]+1j*r.x[d:]; psi=psi/psi[0]
        print('found trial',trial,'gap',r.fun-target,file=sys.stderr)
        with open(f'sic{d}.txt','w') as fo:
            fo.write(' '.join(repr(z.real) for z in psi)+'\n'+' '.join(repr(z.imag) for z in psi)+'\n')
        break
else: print('none',file=sys.stderr)
