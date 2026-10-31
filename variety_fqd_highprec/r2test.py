import mpmath as mp, numpy as np, sys, math
N=int(sys.argv[1]); a=1; mp.mp.dps=int(sys.argv[2]) if len(sys.argv)>2 else 400
z2N=mp.expjpi(mp.mpf(1)/N)
met=lambda x: z2N**(a*(x%N)*((x%N)+N))
half=[u for u in range(1,N) if u<(-u)%N]; selfi=[u for u in range(1,N) if u==(-u)%N]
def full(vals):
    E=lambda x: vals[x%N]
    C=-sum(vals)/mp.sqrt(N)
    return [met(x+y)/(met(x)*met(y))*E(x)*E(y)-C-sum(E(y-z)*E(z)*met(z)*E(x-z) for z in range(N))/mp.sqrt(N) for x in range(N) for y in range(N)]
def build(e,e0,s):
    v=[None]*N; v[0]=e0
    for i,u in enumerate(half): v[u]=e[i]; v[(-u)%N]=1/(met(u)*e[i])
    for j,u in enumerate(selfi): v[u]=s[j]
    return v
def refine(sol):
    e0=mp.mpf(sol[0].real); e0=(mp.sqrt(N)+(1 if e0>0 else -1)*mp.sqrt(N+4))/2
    s=[ (1 if (sol[u]*complex(met(u))**0.5).real>0 else -1)/mp.sqrt(met(u)) for u in selfi]
    # fix sign of sqrt for self-inverse
    s=[ min([si,-si],key=lambda t: abs(complex(t)-sol[u])) for si,u in zip(s,selfi)]
    e=[mp.mpc(sol[u].real,sol[u].imag) for u in half]
    for it in range(60):
        F=mp.matrix(full(build(e,e0,s)))
        if mp.norm(F)<mp.mpf(10)**(-mp.mp.dps+20): break
        h=mp.mpf(10)**(-mp.mp.dps//2); J=mp.matrix(len(F),len(e))
        for k in range(len(e)):
            e2=list(e); e2[k]+=h; F2=mp.matrix(full(build(e2,e0,s)))
            for r in range(len(F)): J[r,k]=(F2[r]-F[r])/h
        JH=J.H; de=mp.lu_solve(JH*J,-(JH*F)); e=[e[k]+de[k] for k in range(len(e))]
    return build(e,e0,s), mp.norm(F)
def inL(v):
    z4N=mp.expjpi(mp.mpf(1)/(2*N)); phi=sum(1 for k in range(1,4*N+1) if math.gcd(k,4*N)==1)
    sq=int(math.isqrt(N+4)); ext=[1] if sq*sq==N+4 else [1,mp.sqrt(N+4)]
    basis=[z4N**j*t for t in ext for j in range(phi)]
    lam=mp.pi
    vec=[mp.re(v)+lam*mp.im(v)]+[mp.re(b)+lam*mp.im(b) for b in basis]
    rel=mp.pslq(vec,maxcoeff=10**40,maxsteps=10**6)
    if rel is None or rel[0]==0: return False, None
    chk=rel[0]*v+sum(r*b for r,b in zip(rel[1:],basis))
    return abs(chk)<mp.mpf(10)**(-mp.mp.dps//2), max(abs(r) for r in rel)
S=[]
for line in open(f'fqd_N{N}_a{a}.txt'):
    p=line.split('|'); S.append(np.array(p[0].split(),float)+1j*np.array(p[1].split(),float))
for k,sol in enumerate(S):
    if sol[0].real<0: continue
    vals,res=refine(sol)
    if res>mp.mpf(10)**(-50): print(k,'refine failed',res); continue
    E1=vals[1]; Em1=vals[N-1]
    a1,h1=inL(E1); a2,h2=inL(E1+Em1)
    print(f'N={N} sol{k}: residual {mp.nstr(res,3)}; E(1) in L? {a1}; E(1)+E(-1) in L? {a2} (height {h2})')
