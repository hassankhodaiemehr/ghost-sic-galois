import mpmath as mp, json, sys, numpy as np
k=sys.argv[1]; mp.mp.dps=int(sys.argv[2])
D=json.load(open(f'scan/g{k}.json')); N=D['N']; neg=D['neg']; add=D['add']
met=[m[0]*mp.expjpi(mp.mpf(m[1])/N) for m in D['met']]
bic=[[mp.expjpi(mp.mpf(2*D['bic'][i][j])/N) for j in range(N)] for i in range(N)]
sub=lambda x,y: add[x][neg[y]]
lines=open(f'sol_{k}.txt').read().split('\n'); zero=int(lines[0])-1
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
            for q in range(N): r[q]+=1/sN
            for z in range(N):
                a,b,c=sub(y,z),z,sub(x,z); w=met[z]/sN
                r[a]-=w*E[b]*E[c]; r[b]-=w*E[a]*E[c]; r[c]-=w*E[a]*E[b]
            rows.append(r)
    return mp.matrix(rows)
def refine(sol):
    E=[mp.mpc(v.real,v.imag) for v in sol]; nr=None
    for it in range(60):
        Fv=mp.matrix(F(E)); nr=mp.norm(Fv)
        if nr<mp.mpf(10)**(-mp.mp.dps+20): break
        Jm=J(E); JH=Jm.H; d=mp.lu_solve(JH*Jm,-(JH*Fv)); E=[E[q]+d[q] for q in range(N)]
    return E,nr
out=open(f't_{k}.txt','w')
reps=[u for u in range(N) if u!=zero and u<=neg[u]]
for i,s in enumerate([s for s in S if s[zero].real>0]):
    E,nr=refine(s)
    for u in reps:
        su=mp.sqrt(met[u]); x=su*E[u]; t=x+1/x
        out.write(f"{k} {N} {i} {u} {mp.nstr(mp.re(t),mp.mp.dps-15)} {mp.nstr(abs(mp.im(t)),3)} {mp.nstr(nr,3)}\n")
out.close(); print(k,N,len(S),'done')
