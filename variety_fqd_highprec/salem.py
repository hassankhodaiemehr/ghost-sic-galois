import mpmath as mp, numpy as np, sys
sys.argv=[sys.argv[0],sys.argv[1],sys.argv[2]]
N=int(sys.argv[1]); mp.mp.dps=int(sys.argv[2])
exec(open('r2test.py').read().split('def inL')[0].replace("N=int(sys.argv[1]); a=1; mp.mp.dps=int(sys.argv[2]) if len(sys.argv)>2 else 400","a=1"))
S=[]
for line in open(f'fqd_N{N}_a1.txt'):
    p=line.split('|'); S.append(np.array(p[0].split(),float)+1j*np.array(p[1].split(),float))
sol=[s for s in S if s[0].real>0][0]
vals,res=refine(sol)
print('N',N,'residual',mp.nstr(res,3))
def sq(u): # consistent sqrt of <u>, even in u
    v=mp.sqrt(met(u)); return v
for u in range(1,(N+1)//2):
    x=sq(u)*vals[u]; t=x+1/x
    poly=None
    for deg in range(1,41):
        p=mp.findpoly(mp.re(t),deg,maxcoeff=10**12,maxsteps=200000) if abs(mp.im(t))<mp.mpf(10)**(-50) else None
        if p: poly=p; break
    if not poly: print(' u',u,'t not found up to deg 40 (Im t=',mp.nstr(mp.im(t),3),')'); continue
    roots=mp.polyroots(poly,maxsteps=200,extraprec=400)
    nr=sum(1 for r in roots if abs(mp.im(r))<mp.mpf(10)**(-30))
    big=sum(1 for r in roots if abs(mp.im(r))<mp.mpf(10)**(-30) and abs(mp.re(r))>2)
    print(f' u={u}: t=x+1/x has deg {len(poly)-1}, real roots {nr}/{len(poly)-1}, with |t|>2: {big}, in [-2,2]: {nr-big}; poly={poly}')
