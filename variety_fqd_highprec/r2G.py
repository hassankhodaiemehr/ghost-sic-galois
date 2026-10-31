import mpmath as mp, json, sys, math, numpy as np, itertools
tag=sys.argv[1]; mp.mp.dps=int(sys.argv[2]); cond=int(sys.argv[3])   # cond: K0=Q(zeta_cond)
D=json.load(open(tag+'.json')); N=D['N']; neg=D['neg']; add=D['add']
met=[m[0]*mp.expjpi(mp.mpf(m[1])/N) for m in D['met']]
bic=[[mp.expjpi(mp.mpf(2*D['bic'][i][j])/N) for j in range(N)] for i in range(N)]
sub=lambda x,y: add[x][neg[y]]
lines=open(f'fqdW_{tag}.txt').read().split('\n'); zero=int(lines[0])-1
S=[np.array(l.split('|')[0].split(),float)+1j*np.array(l.split('|')[1].split(),float) for l in lines[1:] if l.strip()]
sN=mp.sqrt(N)
def F(E):
    C=-sum(E)/sN
    return [bic[x][y]*E[x]*E[y]-C-sum(E[sub(y,z)]*E[z]*met[z]*E[sub(x,z)] for z in range(N))/sN for x in range(N) for y in range(N)]
def J(E):
    rows=[]
    for x in range(N):
        for y in range(N):
            r=[mp.mpc(0)]*N
            r[x]+=bic[x][y]*E[y]; r[y]+=bic[x][y]*E[x]
            for k in range(N): r[k]+=1/sN
            for z in range(N):
                a,b,c=sub(y,z),z,sub(x,z); w=met[z]/sN
                r[a]-=w*E[b]*E[c]; r[b]-=w*E[a]*E[c]; r[c]-=w*E[a]*E[b]
            rows.append(r)
    return mp.matrix(rows)
def refine(sol):
    E=[mp.mpc(v.real,v.imag) for v in sol]
    for it in range(40):
        Fv=mp.matrix(F(E)); nr=mp.norm(Fv)
        if nr<mp.mpf(10)**(-mp.mp.dps+15): break
        Jm=J(E); JH=Jm.H; d=mp.lu_solve(JH*Jm,-(JH*Fv)); E=[E[k]+d[k] for k in range(N)]
    return E,nr
phi=sum(1 for k in range(1,cond+1) if math.gcd(k,cond)==1)
zc=mp.expjpi(mp.mpf(2)/cond)
sq=math.isqrt(N+4); ext=[1] if sq*sq==N+4 else [1,mp.sqrt(N+4)]
basis=[zc**j*t for t in ext for j in range(phi)]
def inL(v):
    vec=[mp.re(v)+mp.pi*mp.im(v)]+[mp.re(b)+mp.pi*mp.im(b) for b in basis]
    rel=mp.pslq(vec,maxcoeff=10**50,maxsteps=10**6)
    if rel is None or rel[0]==0: return False
    return abs(rel[0]*v+sum(r*b for r,b in zip(rel[1:],basis)))<mp.mpf(10)**(-mp.mp.dps//2)
Ref=[]
for s in S:
    if s[zero].real>0:
        E,nr=refine(s); Ref.append(E); print('refined residual',mp.nstr(nr,3),flush=True)
u=[i for i in range(N) if i!=zero and neg[i]!=i][0]
rr=[1,3,-2,5,7,-4,2,11,-3,6,1,-5,9,2,-7,4][:N]; vals=[sum(r*E[x] for r,x in zip(rr,range(N))) for E in Ref]
print('coordinate u =',u,'; deg-1 test (E(u) in L):',[inL(v) for v in vals],flush=True)
for (i,j) in itertools.combinations(range(len(Ref)),2):
    isrefl=all(abs(Ref[j][x]-Ref[i][neg[x]])<1e-20 for x in range(N))
    print(f'pair ({i},{j}) reflection-pair={isrefl}: E_i(u)+E_j(u) in L? {inL(vals[i]+vals[j])}; product in L? {inL(vals[i]*vals[j])}',flush=True)
print('sum of all four in L?',inL(sum(vals)))
