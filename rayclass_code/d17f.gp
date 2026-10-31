read("d17polys.gp");
default(realprecision,20500);
read("g17x.gp");
{
my(e=real(g17[p17]), u=e^2, t=u+1/u, Q=Q17, Dsf=7, dp=17, bnf=bnfinit(y^2-Dsf,1), nf=bnf.nf, rts=nf.roots, idp=if(rts[1]>0,1,2), oth=3-idp, fa, Rt, best=+oo);
my(r=polroots(Q)); print("deg ",poldegree(Q)," nonreal ",#select(w->abs(imag(w))>1e-30,r)," in(-2,2) ",#select(w->abs(imag(w))<1e-30&&abs(real(w))<2,r)," P(0) ",if(polcoef(Q,0)%2,"odd","even"));
fa=nffactor(nf,Q); print("factors over F: ",#fa[,1]);
for(i=1,#fa[,1], my(v=abs(subst(subst(lift(fa[i,1]),'y,rts[idp]),'x,t))); if(v<best,best=v;Rt=lift(fa[i,1])));
my(c=rnfconductor(bnf,Rt)); print("t-field over F: conductor ",idealfactor(nf,c[1][1])," inf ",c[1][2]," Gal ",c[2].cyc," H ",c[3]);
my(Ru=simplify('x^poldegree(Rt,'x)*subst(Rt,'x,'x+1/'x)));
my(bm=bnrinit(bnf,[dp,vector(2,j,j==oth)]), cyc=bm.cyc, ord(v)=lcm(vector(#v,i,cyc[i]/gcd(v[i],cyc[i]))));
print("Cl(17 inf-) = ",cyc,"   Cl(17) = ",bnrinit(bnf,dp).cyc);
my(beta, found=0); for(a=-200,200, for(b=-200,200, if(found,break); my(g=1+dp*(a+b*'y), s1=subst(g,'y,rts[idp]), s2=subst(g,'y,rts[oth])); if(s1>0&&s2<0, beta=g; found=1)));
my(cb=bnrisprincipal(bm,idealhnf(nf,beta),0)); print("class of c_-: ",cb," order ",ord(cb));
my(nt=0, ok=0, bad=0, kap_ok=0, kap_n=0);
forprime(p=3,30000, if(p==2||p==3||p==7||p==17, next);
  my(dec=idealprimedec(nf,p));
  for(k=1,#dec, my(P=dec[k]); if(P.f!=1, next);
    my(modP=nfmodprinit(nf,P), Rp=Pol(apply(c->nfmodpr(nf,c,modP),Vec(Ru)),'x));
    if(poldisc(Rp)==0, next);
    my(fr=factor(Rp), degs=apply(poldegree,fr[,1]), cl=bnrisprincipal(bm,P,0), o=ord(cl));
    nt++; if(vecmin(degs)==o && vecmax(degs)==o, ok++, bad++);
    if(prod(i=1,#cl,(cl[i]-cb[i])%cyc[i]==0),
       kap_n++; if(vecmax(degs)==2 && vecmin(degs)==2 && #select(f->polcoef(f,0)!=1,fr[,1])==0, kap_ok++))));
print("degree-1 primes tested: ",nt,";  Frobenius order = ray-class order: ",ok,"  mismatches: ",bad);
print("primes in class of c_-: ",kap_n,";  Frobenius acts as kappa: ",kap_ok);
my(a=subst(Q,'x,2), b=subst(Q,'x,-2)); print("P(2) = ",factor(a,10^8)); print("P(-2) = ",factor(b,10^8));
}
