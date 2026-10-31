"""Reduced finite-quantum-dilogarithm system for the SIC metric groups G_d.

Conventions verified against Radchenko-Wheeler (arXiv:2609.21892v2), Thm. 5:
  alpha = [[0,-1],[1,d-1]]  (tr = d-1),  gamma = alpha^3
  gamma - I = d*alpha*(alpha-I) = d((d-2)alpha - I)      [exact, all d]
  G_d = Z^2 / Z^2(gamma-I)  ~=  Z/d x Z/(d(d-3)),  N = |G_d| = d^2(d-3)
  <u>   = (-1)^(u1+u2+u1u2) * e_{2N}(omega(u gamma, u)),  omega(v,w)=v1w2-v2w1
  <x;y> = e_N(omega(x gamma - x, y)) = <x+y>/(<x><y>)
  F has kernel <x;y>/sqrt(N);  F^4 = I,  F^2 = rho  (rho f(u)=f(-u))
  W = <.>^{-1} F^{-1};  W^3 = c*I with |c| = 1
  lambda^3 = N^{-1/2} sum_x <x>^{-1};  system: W E = lambda E,
             E(u)E(-u) = <u>^{-1} (u != 0),  E(0)^2 = sqrt(N) E(0) + 1

In x-coordinates x_u = s_u E(u), s_u^2 = <u>, s_{-u}=s_u, the pairing is x_u x_{-u} = 1.
Splitting by rho and writing t_u = x_u + x_{-u}, y_u = x_u - x_{-u} gives the
equivalent form used in the analysis:

      t in 2 X_lambda^+ ,   y in X_lambda^- ,   y_u^2 = t_u^2 - 4 .

All reality content sits in the last (scalar, one per pair) equation.
"""
from math import gcd
import numpy as np

def xgcd(a,b):
    old_r,r=a,b; old_s,s=1,0; old_t,t=0,1
    while r: q=old_r//r; old_r,r=r,old_r-q*r; old_s,s=s,old_s-q*s; old_t,t=t,old_t-q*t
    if old_r<0: old_r,old_s,old_t=-old_r,-old_s,-old_t
    return old_r,old_s,old_t

def matmul(X,Y): return [[sum(X[i][k]*Y[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
def omega(v,w): return v[0]*w[1]-v[1]*w[0]
def rowmat(u,M): return (u[0]*M[0][0]+u[1]*M[1][0], u[0]*M[0][1]+u[1]*M[1][1])

def row_hnf(M):
    """Upper-triangular basis of the lattice spanned by the ROWS of 2x2 M."""
    r1,r2=list(M[0]),list(M[1])
    g,a,b=xgcd(r1[0],r2[0])
    A=[a*r1[0]+b*r2[0], a*r1[1]+b*r2[1]]
    B=[(r2[0]//g)*r1[0]-(r1[0]//g)*r2[0], (r2[0]//g)*r1[1]-(r1[0]//g)*r2[1]]
    assert B[0]==0
    h22=abs(B[1])
    return [[g, A[1]%h22],[0,h22]]

def metric_group(d):
    A=[[0,-1],[1,d-1]]; G3=matmul(matmul(A,A),A)
    M=[[G3[0][0]-1,G3[0][1]],[G3[1][0],G3[1][1]-1]]
    N=abs(M[0][0]*M[1][1]-M[0][1]*M[1][0]); assert N==d*d*(d-3)
    H=row_hnf(M); h11,h12,h22=H[0][0],H[0][1],H[1][1]; assert h11*h22==N
    def red(u):
        k=u[0]//h11; u=(u[0]-k*h11,u[1]-k*h12); return (u[0],u[1]%h22)
    reps=[(i,j) for i in range(h11) for j in range(h22)]
    return dict(d=d,N=N,gamma=G3,M=M,reps=reps,idx={r:k for k,r in enumerate(reps)},red=red)

def operators(d):
    B=metric_group(d); N=B['N']; reps=B['reps']; G3=B['gamma']; red=B['red']; idx=B['idx']
    q=lambda u:(-1)**(u[0]+u[1]+u[0]*u[1])*np.exp(2j*np.pi*omega(rowmat(u,G3),u)/(2*N))
    bc=lambda x,y: np.exp(2j*np.pi*omega((rowmat(x,G3)[0]-x[0],rowmat(x,G3)[1]-x[1]),y)/N)
    F=np.array([[bc(x,y) for y in reps] for x in reps])/np.sqrt(N)
    Dg=np.diag([q(u) for u in reps])
    P=np.zeros((N,N))
    for a,u in enumerate(reps): P[a, idx[red((-u[0],-u[1]))]]=1
    W=np.linalg.inv(Dg)@np.linalg.inv(F)
    lam3=sum(1/q(x) for x in reps)/np.sqrt(N)
    B.update(q=q,bc=bc,F=F,D=Dg,rho=P,W=W,lam3=lam3)
    return B

def dimensions(d):
    B=operators(d); N=B['N']; W=B['W']; P=B['rho']; lam3=B['lam3']
    ev=np.linalg.eigvals(W)
    cub=[abs(lam3)**(1/3)*np.exp(1j*(np.angle(lam3)+2*np.pi*k)/3) for k in range(3)]
    lam=min(cub,key=lambda c: np.abs(ev-c).min())
    dimX=int((np.abs(ev-lam)<1e-6).sum())
    _,s,Vh=np.linalg.svd(W-lam*np.eye(N)); K=Vh[N-dimX:].conj().T
    w=np.linalg.eigvals(np.linalg.lstsq(K,P@K,rcond=None)[0])
    tors=[u for u in B['reps'] if B['red']((2*u[0],2*u[1]))==(0,0) and u!=(0,0)]
    return dict(d=d,N=N,dimX=dimX,t_dim=int((np.abs(w-1)<1e-6).sum()),
                y_dim=int((np.abs(w+1)<1e-6).sum()),
                pairs=(N-1-len(tors))//2,two_torsion=len(tors))

if __name__=='__main__':
    for d in range(4,10): print(dimensions(d))
