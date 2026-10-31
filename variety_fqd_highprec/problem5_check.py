import numpy as np
from scipy.optimize import minimize
d=4
def Wgamma(alpha):
    # partial transpose of (I + alpha*V), V^Gamma = d * |phi+><phi+|
    phi=np.eye(d).reshape(d*d)/np.sqrt(d)
    return np.eye(d*d)+alpha*d*np.outer(phi,phi)
def two_copy(alpha):
    G=Wgamma(alpha).reshape(d,d,d,d)  # a b a' b'
    # order A1 A2 B1 B2
    T=np.einsum('abcd,efgh->aebfcgdh',G,G).reshape(d*d*d*d,d*d*d*d)
    return T
def run(alpha,trials=40,rng=np.random.default_rng(0)):
    T=two_copy(alpha); D=d*d
    def f(x):
        z=x[:4*D]+1j*x[4*D:]
        U=z[:2*D].reshape(2,D); V=z[2*D:].reshape(2,D)
        psi=np.einsum('ia,ib->ab',U,V).reshape(-1)
        n=np.vdot(psi,psi).real
        return (np.vdot(psi,T@psi).real)/n
    best=np.inf
    for _ in range(trials):
        r=minimize(f,rng.normal(size=8*D),method='L-BFGS-B')
        best=min(best,r.fun)
    return best
for a in [-0.45,-0.5,-0.52,-0.55]:
    print(a, run(a))
