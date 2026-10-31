import mpmath as mp, numpy as np, sys
d=int(sys.argv[1]); mp.mp.dps=int(sys.argv[2])
P=[]
for line in open(f'sols_d{d}.txt'):
    re_,im_=line.split('|'); P.append(np.array(re_.split(),float)+1j*np.array(im_.split(),float))
w=np.exp(2j*np.pi/d); X=np.roll(np.eye(d),1,axis=0); Zm=np.diag(w**np.arange(d)); UP=np.eye(d)[[(-j)%d for j in range(d)]]
same=lambda u,v: abs(abs(np.vdot(u/np.linalg.norm(u),v/np.linalg.norm(v)))-1)<1e-7
gh=[v for v in P if not same(np.conj(v[d:]),v[:d]) and same(np.conj(v[d:]),UP@v[:d])]
print('U_P ghosts:',len(gh),file=sys.stderr)
v=gh[int(sys.argv[3]) if len(sys.argv)>3 else 0]
# mp matrices
om=mp.expjpi(mp.mpf(2)/d); tau=-mp.expjpi(mp.mpf(1)/d)
Xm=mp.matrix(d,d); Zmm=mp.matrix(d,d)
for j in range(d): Xm[(j+1)%d,j]=1; Zmm[j,j]=om**j
def D(a,b): return (tau**(a*b))*(Xm**a)*(Zmm**b)
ps=[(a,b) for a in range(d) for b in range(d) if (a,b)!=(0,0)]
Ds={p:D(*p) for p in ps}; DH={p:Ds[p].H for p in ps}
x=[mp.mpc(c) for c in v[:d]]; y=[mp.mpc(c) for c in v[d:]]
x=[c/x[0] for c in x]; s=sum(a*b for a,b in zip(y,x)); y=[c/s for c in y]
def bil(M,y,x): return sum(y[i]*M[i,j]*x[j] for i in range(d) for j in range(d))
def F(z):
    x=z[:d]; y=z[d:]; s=sum(a*b for a,b in zip(y,x))
    return [(d+1)*bil(Ds[p],y,x)*bil(DH[p],y,x)-s**2 for p in ps]+[x[0]-1, s-1]
def J(z,h=None):
    # analytic Jacobian
    x=z[:d]; y=z[d:]; s=sum(a*b for a,b in zip(y,x)); rows=[]
    for p in ps:
        A=Ds[p]; B=DH[p]; a=bil(A,y,x); b=bil(B,y,x)
        r=[]
        for k in range(d):
            da=sum(y[i]*A[i,k] for i in range(d)); db=sum(y[i]*B[i,k] for i in range(d))
            r.append((d+1)*(da*b+a*db)-2*s*y[k])
        for k in range(d):
            da=sum(A[k,j]*x[j] for j in range(d)); db=sum(B[k,j]*x[j] for j in range(d))
            r.append((d+1)*(da*b+a*db)-2*s*x[k])
        rows.append(r)
    rows.append([1]+[0]*(2*d-1)); rows.append([y[k] for k in range(d)]+[x[k] for k in range(d)])
    return mp.matrix(rows)
z=x+y
for it in range(50):
    Fv=mp.matrix(F(z)); nr=mp.norm(Fv)
    if nr<mp.mpf(10)**(-mp.mp.dps+15): break
    Jm=J(z); dz=mp.lu_solve(Jm.H*Jm,-(Jm.H*Fv)); z=[z[k]+dz[k] for k in range(2*d)]
print('residual',mp.nstr(nr,3),file=sys.stderr)
x=z[:d]; y=z[d:]; s=sum(a*b for a,b in zip(y,x))
for p in ps:
    e=mp.sqrt(d+1)*bil(Ds[p],y,x)/s
    print(p[0],p[1],mp.nstr(mp.re(e),mp.mp.dps-20),mp.nstr(mp.im(e),5))
