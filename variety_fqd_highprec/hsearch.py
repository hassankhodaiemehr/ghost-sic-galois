import numpy as np, itertools, sys, time
from scipy.optimize import minimize

class Heis:
    """Weyl-Heisenberg group over finite abelian group A = Z_n1 x ... x Z_nk."""
    def __init__(self, shape):
        self.shape=tuple(shape); self.d=int(np.prod(shape))
        self.shifts=list(itertools.product(*[range(n) for n in shape]))
    def overlaps(self, psi):
        P=psi.reshape(self.shape); d=self.d
        C=np.empty((d,)+self.shape, complex)
        for i,a in enumerate(self.shifts):
            u=np.conj(np.roll(P,[-x for x in a],axis=tuple(range(len(a)))))*P   # conj(psi(y+a)) psi(y)
            C[i]=d*np.fft.ifftn(u)
        return C   # C[a][b] = <psi|X^a Z^b|psi> up to phase
    def f_and_grad(self, x):
        d=self.d; psi=x[:d]+1j*x[d:]
        n=np.vdot(psi,psi).real
        C=self.overlaps(psi); A2=np.abs(C)**2
        S=np.sum(A2**2)
        f=S/n**4
        P=psi.reshape(self.shape); M=np.zeros(self.shape,complex)
        for i,a in enumerate(self.shifts):
            w=A2[i]*np.conj(C[i])
            W=d*np.fft.ifftn(w)                    # W_a(y) = sum_b w_ab chi_b(y)
            M+=np.roll(W*P,list(a),axis=tuple(range(len(a))))  # wait: sign conventions checked numerically
        g=4*M.reshape(-1)/n**4 - 4*S/n**5*psi
        return f, np.concatenate([2*g.real,2*g.imag])
    def target(self):
        d=self.d; return 2*d/(d+1)

def search(shape, starts=50, seed=0, tol=1e-10, verbose=False):
    H=Heis(shape); d=H.d; rng=np.random.default_rng(seed)
    best=(np.inf,None); hits=0
    for s in range(starts):
        x0=rng.normal(size=2*d)
        r=minimize(H.f_and_grad,x0,jac=True,method='L-BFGS-B',options={'maxiter':20000,'ftol':1e-16,'gtol':1e-12})
        gap=r.fun-H.target()
        if gap<tol: hits+=1
        if r.fun<best[0]: best=(r.fun,r.x)
    return H, best[0]-H.target(), hits, best[1]

if __name__=='__main__':
    # gradient check
    H=Heis((2,3)); x=np.random.default_rng(1).normal(size=12)
    f,g=H.f_and_grad(x); e=1e-6; num=np.array([(H.f_and_grad(x+e*np.eye(12)[i])[0]-H.f_and_grad(x-e*np.eye(12)[i])[0])/(2*e) for i in range(12)])
    print('grad check', np.max(np.abs(num-g)))
