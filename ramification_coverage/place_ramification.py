
"""Real-quadratic arithmetic for the SIC / place-law program.

d        SIC dimension (>=4)
D        = (d-3)(d+1) = (d-1)^2 - 4
D0       squarefree part of D;  F = Q(sqrt(D0))
O_F      = Z[w],  w = (1+sqrt D0)/2 if D0 = 1 (mod 4), else w = sqrt(D0)
dprime   = d (d odd), 2d (d even)
eps_a    = (d-1 + sqrt D)/2   (eigenvalue of alpha, tr alpha = d-1, N(eps_a)=+1)
"""
from math import isqrt

def factorize(n):
    f={}
    for p in (2,3,5):
        while n%p==0: f[p]=f.get(p,0)+1; n//=p
    p,step,i=7,[4,2,4,2,4,6,2,6],0
    while p*p<=n:
        while n%p==0: f[p]=f.get(p,0)+1; n//=p
        p+=step[i]; i=(i+1)%8
    if n>1: f[n]=f.get(n,0)+1
    return f

def squarefree_part(n):
    s=1
    for p,e in factorize(n).items():
        if e%2: s*=p
    return s

def is_squarefree(n):
    return all(e==1 for e in factorize(n).values())

def icbrt(n):
    if n<0: return -icbrt(-n)
    r=int(round(n**(1/3))) if n < 10**15 else int(round(n**(1/3)))
    while r**3>n: r-=1
    while (r+1)**3<=n: r+=1
    return r

# ---------- Pell: minimal x^2 - D0 y^2 = +-1 via CF of sqrt(D0) ----------
def pell_pm1(D0):
    a0=isqrt(D0); m,dd,a=0,1,a0
    num1,num=1,a0; den1,den=0,1
    while True:
        if num*num - D0*den*den in (1,-1):
            return num,den,num*num-D0*den*den
        m=dd*a-m; dd=(D0-m*m)//dd; a=(a0+m)//dd
        num1,num=num,a*num+num1
        den1,den=den,a*den+den1

# ---------- fundamental unit of O_F as (a,b):  eps=(a+b sqrt D0)/2 ----------
def fundamental_unit_ab(D0):
    x,y,nrm=pell_pm1(D0)
    a,b=2*x,2*y                      # eps = x + y sqrt D0 = (2x+2y sqrt)/2
    if D0%4==1:
        # O_F may contain a cube root of eps:  a(a^2 -+ 3) = 2x
        for s,target in ((1,2*x),(-1,2*x)):
            pass
        for sign in (1,-1):          # sign = N(eps_small)
            # N=+1 : a(a^2-3)=2x ;  N=-1 : a(a^2+3)=2x
            c = -3 if sign==1 else 3
            r=icbrt(2*x)
            for aa in range(max(1,r-2),r+3):
                if aa*(aa*aa+c)==2*x:
                    num=aa*aa-4*sign
                    if num>0 and num%D0==0:
                        bb=isqrt(num//D0)
                        if bb*bb*D0==num and (aa-bb)%2==0 and aa*aa-D0*bb*bb==4*sign:
                            return aa,bb,sign
    return a,b,nrm

# ---------- Z[w] arithmetic ----------
class Fld:
    def __init__(self,D0):
        self.D0=D0
        if D0%4==1: self.Tw,self.Nw,self.half,self.disc=1,(1-D0)//4,True,D0
        else:       self.Tw,self.Nw,self.half,self.disc=0,-D0,False,4*D0
    def mul(self,x,y,n=None):
        x0,x1=x; y0,y1=y
        c0=x0*y0-x1*y1*self.Nw
        c1=x0*y1+x1*y0+x1*y1*self.Tw
        return (c0%n,c1%n) if n else (c0,c1)
    def pw(self,x,e,n):
        r=(1%n,0); b=(x[0]%n,x[1]%n)
        while e:
            if e&1: r=self.mul(r,b,n)
            b=self.mul(b,b,n); e>>=1
        return r
    def norm(self,x):
        x0,x1=x
        return x0*x0+x0*x1*self.Tw+x1*x1*self.Nw
    def to_w(self,a,b):
        """(a+b sqrt D0)/2  ->  (c0,c1) with value c0+c1 w."""
        if self.half:
            assert (a-b)%2==0
            return ((a-b)//2, b)
        assert a%2==0 and b%2==0
        return (a//2, b//2)

def splitting(disc,p):
    if p==2:
        if disc%2: return 1 if disc%8==1 else -1
        return 0
    if disc%p==0: return 0
    return 1 if pow(disc%p,(p-1)//2,p)==1 else -1

def unit_group_order(F,n):
    tot=1
    for p,e in factorize(n).items():
        s=splitting(F.disc,p)
        if s==1:    tot*= p**(2*e-2)*(p-1)*(p-1)
        elif s==-1: tot*= p**(2*e-2)*(p*p-1)
        else:       tot*= p**(2*e-1)*(p-1)
    return tot

def mult_order(F,x,n,bound):
    r=bound
    for p in factorize(bound):
        while r%p==0 and F.pw(x,r//p,n)==(1%n,0):
            r//=p
    assert F.pw(x,r,n)==(1%n,0)
    return r


# ======================================================================
# Driver: verification of (C_d), of ord(eps_alpha), and of the
# Richaud-Degert family, for 4 <= d <= DMAX.
#
#   (C_d)  no unit of F congruent to 1 mod d' is negative at infinity_-
#          <=>  infinity_- ramifies in the ray class field of conductor d' inf_-
#          <=>  |Cl(d' inf_-)| = 2 |Cl(d')|
# ======================================================================
from math import isqrt
import csv, sys

def analyse(d):
    D  = (d-3)*(d+1)
    D0 = squarefree_part(D)
    g  = isqrt(D//D0); assert g*g*D0 == D
    F  = Fld(D0)
    dpr = d if d % 2 else 2*d
    ea = F.to_w(d-1, g)                      # eps_alpha = ((d-1)+g sqrt D0)/2
    assert F.norm(ea) == 1
    a, b, nu = fundamental_unit_ab(D0)
    ef = F.to_w(a, b)                        # fundamental unit, > 1 at infinity_+
    G  = unit_group_order(F, dpr)
    r  = mult_order(F, ef, dpr, G)           # order of eps mod d'
    m1 = ((-1) % dpr, 0)
    has_m1 = (r % 2 == 0) and F.pw(ef, r//2, dpr) == m1   # is -1 in <eps mod d'> ?
    # eps > 0 at infinity_+ ; sgn_{infinity_-}(eps) = sign of N(eps)
    #   (C_d) fails  <=>  some  +- eps^n  is  = 1 mod d'  and  < 0 at infinity_-
    fail_a = (nu == -1 and r % 2 == 1)                       # eta = eps^n, n odd
    fail_b = has_m1 and (nu == 1 or (r//2) % 2 == 0)         # eta = -eps^n
    ra = mult_order(F, ea, dpr, G)                           # order of eps_alpha mod d'
    return dict(d=d, D=D, D0=D0, dpr=dpr, nu=nu, r=r, ord_ea=ra,
                fund=(ea == ef), ram=not (fail_a or fail_b))

def main(DMAX=20000, out='place_ramification_d_le_20000.csv'):
    rows, notram, nfund = [], [], 0
    for d in range(4, DMAX+1):
        x = analyse(d); rows.append(x)
        if not x['ram']: notram.append(d)
        nfund += x['fund']
        # Lemma: eps_alpha^3 = 1 (d odd) / 1-d (d even), modulo d'
        F = Fld(x['D0']); g = isqrt(x['D']//x['D0']); ea = F.to_w(d-1, g)
        c3 = F.pw(ea, 3, x['dpr'])
        assert c3 == ((1 if d % 2 else 1-d) % x['dpr'], 0), d
        assert x['ord_ea'] == (3 if d % 2 else 6), d
    fam = [d for d in range(5, DMAX+1, 2)
           if ((d-1)//2) % 2 == 0 and is_squarefree(((d-1)//2)**2 - 1)]
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['d','D','D0','dpr','nu','r','ord_ea','fund','ram'])
        w.writeheader()
        for x in rows: w.writerow(x)
    n = len(rows)
    print(f"dimensions examined                     : {n}")
    print(f"(C_d) holds (inf_- ramifies)            : {n-len(notram)}   exceptions: {notram}")
    print(f"O_F^x = <-1, eps_alpha>                 : {nfund}  ({nfund/n:.3%})")
    print(f"Richaud-Degert family D (unconditional) : {len(fam)}  ({len(fam)/n:.3%})")
    print(f"  all of D satisfies O_F^x=<-1,eps_a>   : {all(analyse(d)['fund'] for d in fam[:400])}")
    print(f"  first members                         : {fam[:21]}")

if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20000)
