import numpy as np, sys
def load(d):
    P=[]
    for line in open(f'sols_d{d}.txt'):
        re,im=line.split('|'); v=np.array(re.split(),float)+1j*np.array(im.split(),float); P.append(v)
    return np.array(P)
for d in (4,5):
    P=load(d); w=np.exp(2j*np.pi/d)
    X=np.roll(np.eye(d),1,axis=0); Z=np.diag(w**np.arange(d)); UP=np.eye(d)[[(-j)%d for j in range(d)]]
    Ds={(a,b):np.linalg.matrix_power(X,a)@np.linalg.matrix_power(Z,b) for a in range(d) for b in range(d)}
    same=lambda u,v: abs(abs(np.vdot(u/np.linalg.norm(u),v/np.linalg.norm(v)))-1)<1e-7
    nS=nG=nGpar=nGany=0
    for v in P:
        x,y=v[:d],v[d:]
        if same(np.conj(y),x): nS+=1; continue
        nG+=1
        if same(np.conj(y),UP@x): nGpar+=1
        if any(same(np.conj(y),Ds[q]@UP@x) for q in Ds): nGany+=1
    print(f'd={d}: SIC {nS}, non-real {nG}, parity-Hermitian (U_P) {nGpar}, parity-Hermitian for some D_q U_P {nGany}')
