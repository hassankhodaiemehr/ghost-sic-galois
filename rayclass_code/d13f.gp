default(realprecision,20500);
read("d13polys.gp"); read("g13x.gp");
{
my(e=real(g13[14]), u=e^2, t=u+1/u, Q=Q13, bnf=bnfinit(y^2-35,1), nf=bnf.nf, rts=nf.roots, idp=if(rts[1]>0,1,2), oth=3-idp, fa, Rt, best=+oo, t0=getabstime());
fa=nffactor(nf,Q);
print("nffactor done in ",(getabstime()-t0)\1000,"s, factors: ",#fa[,1]);
for(i=1,#fa[,1], my(v=abs(subst(subst(lift(fa[i,1]),'y,rts[idp]),'x,t))); if(v<best,best=v;Rt=lift(fa[i,1])));
write("d13Rt.gp","Rt13=",Rt,";");
my(Ru=simplify('x^48*subst(Rt,'x,'x+1/'x)));
my(bm=bnrinit(bnf,[13,vector(2,j,j==oth)]), cyc=bm.cyc, ord(v)=lcm(vector(#v,i,cyc[i]/gcd(v[i],cyc[i]))));
\\ class of c_-: beta = 1 mod 13, negative at oth, positive at idp
my(beta, found=0); for(a=-80,80, for(b=-80,80, if(found,break); my(g=1+13*(a+b*'y), s1=subst(g,'y,rts[idp]), s2=subst(g,'y,rts[oth])); if(s1>0&&s2<0, beta=g; found=1)));
my(cb=bnrisprincipal(bm,idealhnf(nf,beta),0)); print("class of c_- in Cl(13 inf-): ",cb," order ",ord(cb));
my(nt=0, ok=0, bad=0, kap_ok=0, kap_n=0);
forprime(p=3,20000, if(p==5||p==7||p==13, next);
  my(dec=idealprimedec(nf,p));
  for(k=1,#dec, my(P=dec[k]); if(P.f!=1, next);
    my(modP=nfmodprinit(nf,P), Rp=Pol(apply(c->nfmodpr(nf,c,modP),Vec(Ru)),'x));
    if(poldisc(Rp)==0, next);
    my(fr=factor(Rp), degs=apply(poldegree,fr[,1]), cl=bnrisprincipal(bm,P,0), o=ord(cl));
    nt++; if(vecmin(degs)==o && vecmax(degs)==o, ok++, bad++);
    if((cl-cb)%1==0 && apply(x->x,(cl-cb))==vector(#cl,i,0)~ || prod(i=1,#cl,(cl[i]-cb[i])%cyc[i]==0),
       kap_n++; if(vecmax(degs)==2 && vecmin(degs)==2 && #select(f->polcoef(f,0)!=1,fr[,1])==0, kap_ok++))));
print("degree-1 primes of F tested: ",nt,";  Frobenius order = ray-class order: ",ok,"  mismatches: ",bad);
print("primes in the class of c_-: ",kap_n,";  Frobenius acts as kappa (x -> 1/x): ",kap_ok);
}
