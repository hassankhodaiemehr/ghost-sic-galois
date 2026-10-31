\\ Reciprocity-law test for the squared ghost overlaps.  Inputs: d, DIG, QMAX set before read.
default(realprecision, DIG+50);
{
W=read(Str("eps",d,".txt")); P=W[1]; V=apply(z->z^2, W[2]); np=#P;
Dl=(d-3)*(d+1); D0=core(Dl); k=sqrtint(Dl/D0);
eps=10.^(-DIG/2);
}
idx(p)={my(a=p[1]%d,b=p[2]%d); for(i=1,np, if(P[i][1]==a&&P[i][2]==b, return(i))); 0};
expr(v)=for(j=1,#Cand, if(abs(Cand[j][2]-v)<eps, return(Cand[j][1]))); my(Q=algdep(v,n), fa=factor(Q)[,1], f=0, rr, best=0); for(i=1,#fa, if(abs(subst(fa[i],'x,v))<eps, f=fa[i])); if(f==0, error("algdep failed")); rr=nfroots(Kw, f); for(i=1,#rr, my(z=subst(lift(rr[i]),'w,e0)); if(abs(z-v)<eps, best=subst(lift(rr[i]),'w,'x))); if(best==0, error("no matching root in K")); best;
emb(a,s)= subst(lift(a),y,s*sqrt(D0));
phi(al,Lam)=my(c=Vec(lift(al)), a0, b0, X=(2*Lam^-1-(d-1))/k); if(#c==1, c=[0,c[1]]); b0=c[1]; a0=c[2]; Mod(a0,d)*matid(2)+Mod(b0,d)*X;
{
\\ 1. stabilizer of the e_p in GL2(Z/d)
S=List();
forvec(m=vector(4,i,[0,d-1]), my(G=Mod([m[1],m[2];m[3],m[4]],d),ok=1);
  if(matdet(G)==0,next);
  for(i=1,np, my(q=lift(G*P[i]~)); if(abs(V[idx(q)]-V[i])>eps, ok=0;break));
  if(ok, listput(S,G)));
S=Vec(S); print("d=",d,"  D0=",D0,"  k=",k,"  stabilizer size ",#S);
Lams=[G | G<-S, G^3==1 && G!=1]; print("order-3 elements: ",apply(lift,Lams), " traces ",apply(G->lift(trace(G)),Lams));
}
{
\\ 2. minimal polynomial of e0 and exact expressions
i0=idx([0,1]); e0=V[i0];
my(n0=2*NDEG); Pm=algdep(e0,n0); my(fa=factor(Pm)[,1]); for(i=1,#fa, if(abs(subst(fa[i],'x,e0))<eps, Pm=fa[i]));
n=poldegree(Pm); print("minpoly degree ",n,"  height ~10^",round(log(vecmax(apply(abs,Vec(Pm))))/log(10)));
pr0=polredbest([Pm,10^6],1); my(pr=pr0); Pm0=Pm; Pm=pr[1]; th=subst(lift(modreverse(pr[2])),'x,e0); if(abs(subst(Pm,'x,th))>eps*vecmax(apply(abs,Vec(Pm)))*max(1,abs(th))^poldegree(Pm), error("theta")); e0=th;
print("polredbest: new defining polynomial height ~10^",round(log(vecmax(apply(abs,Vec(Pm))))/log(10)));
K=nfinit([Pm,10^6]); Kw=nfinit([subst(Pm,'x,'w),10^6]); gal=iferr(alarm(2400, nfgaloisconj(K)), E, print("nfgaloisconj: ",E); []); print("#automorphisms ",#gal, "  (Galois iff = degree)");
e0ex=lift(pr0[2]); Cand=vector(#gal,j, my(z=lift(subst(e0ex,'x,Mod(gal[j],Pm)))); [z, subst(z,'x,e0)]); print("candidates from automorphisms: ",#Cand);
Ex=vector(np); 
for(i=1,np, if(Ex[i]!=0, next);
  my(j0=0); for(j=1,i-1, if(abs(V[j]-V[i])<eps, j0=j; break)); if(j0, Ex[i]=Ex[j0]; next);
  Ex[i]=Mod(expr(V[i]),Pm));
print("expressions done");
for(i=1,np, if(abs(subst(lift(Ex[i]),'x,e0)-V[i])>eps, error("bad expr ",i)));
write(Str("expr",d,".gp"), [Pm, Ex, e0]); print("saved expressions");
\\ exactness checks independent of numerics:  e_p e_{-p} = 1, Pm(e_p)=0
for(i=1,np, if(Ex[i]*Ex[idx(-P[i])]!=1, print("WARNING: e_p e_-p != 1 at ",P[i])));
sq=Mod(expr(sqrt(D0)),Pm); if(sq^2!=D0, error("sqrt D0 not exact"));
print("exact: all e_p in Q(e0), e_p e_{-p}=1, sqrt(D0) in Q(e0)");
\\ abelian over F?  check commutation of automorphisms
ab=1; for(i=1,#gal, for(j=1,#gal, if(subst(gal[i],'x,Mod(gal[j],Pm))!=subst(gal[j],'x,Mod(gal[i],Pm)), ab=0))); if(#gal, print("Gal(K/Q) abelian: ",ab));
fixF=[g | g<-gal, subst(lift(sq),'x,Mod(g,Pm))==sq]; print("|Gal(K/F)| = ",#fixF);
abF=1; for(i=1,#fixF, for(j=1,#fixF, if(subst(fixF[i],'x,Mod(fixF[j],Pm))!=subst(fixF[j],'x,Mod(fixF[i],Pm)), abF=0))); print("Gal(K/F) abelian: ",abF);
kap=[g | g<-fixF, subst(lift(Ex[i0]),'x,Mod(g,Pm))==1/Ex[i0]]; print("kappa (e0->1/e0) in Gal(K/F): ",#kap);
if(#kap, my(g=kap[1]); okk=1; for(i=1,np, if(subst(lift(Ex[i]),'x,Mod(g,Pm))!=1/Ex[i], okk=0)); print("  kappa inverts every e_p: ",okk));
}
{
\\ 3. Frobenius tests
bnf=bnfinit(y^2-D0,1); w=if(D0%4==1,(1+y)/2,y);
fu=bnf.fu[1]; print("fundamental unit ",fu,"  norm ",norm(Mod(fu,y^2-D0)));
units=[Mod(1,y^2-D0),Mod(-1,y^2-D0),Mod(fu,y^2-D0),Mod(-fu,y^2-D0)];
stats=Map(); tot=0; ncoll=0;
den=lcm(vector(np,i,denominator(content(lift(Ex[i])))));den=lcm(den,denominator(content(lift(sq)))); print("expression denominators: ",factor(den));
forprime(q=3,QMAX, if(d%q==0||Dl%q==0||poldisc(Pm)%q==0||den%q==0, next);
  my(dec=idealprimedec(bnf,q));
  foreach(dec,pr,
    my(Nq=pr.p^pr.f, gen=bnfisprincipal(bnf,pr,1)[2], al=Mod(nfbasistoalg(bnf,gen),y^2-D0));
    \\ residue root of sqrt(D0) mod pr (if split) / choose factor g of Pm mod q above pr
    my(fa=factormod(Pm,q)[,1], g=0, perm=vector(np));
    for(t=1,#fa, my(gg=lift(fa[t]), r=Mod(Mod(1,q)*lift(sq),Mod(1,q)*gg));
       \\ sqrt D0 in residue field; pr above if (sqrt D0 - root) in pr : test via valuation of al? use idealval on element sqrt - lift
       if(pr.f==2, g=gg; break);
       my(rv=lift(lift(r))); if(poldegree(rv)>0, next); rv=polcoef(rv,0);
       if(idealval(bnf, y-rv, pr)>0, g=gg; break));
    if(g==0, print("no factor for ",q); next);
    my(red=vector(np,j,Mod(Mod(1,q)*lift(Ex[j]),Mod(1,q)*g)));
    my(coll=0); for(i=1,np, for(j=i+1,np, if(Ex[i]!=Ex[j] && red[i]==red[j], coll=1))); if(coll, ncoll++; next);
    for(i=1,np, my(v=red[i]^Nq); for(j=1,np, if(red[j]==v, perm[i]=j; break)));
    if(vecmin(perm)==0, print("perm incomplete q=",q); next);
    tot++;
    foreach([1,2],li, foreach([1,-1],s,
      my(Lam=if(li==1,Lams[1],Lams[1]^-1), okany=0, okall=1, nchoices=0);
      foreach(units,u, my(pi=al*u); if(emb(pi,s)>0, nchoices++; my(G=phi(pi,Lam), ok=1);
          for(i=1,np, if(perm[i]!=idx(lift(G*P[i]~)) && Ex[perm[i]]!=Ex[idx(lift(G*P[i]~))], ok=0; break));
          okany+=ok; okall*=ok));
      my(key=Str(if(pr.f==1,"split","inert")," Lam",if(li==1,"","^-1")," pos@",if(s==1,"inf+","inf-")));
      mapput(stats,key,if(mapisdefined(stats,key),mapget(stats,key),[0,0,0])+[1,okall,okany])))));
print("prime ideals tested: ",tot,"   (skipped for residue collisions: ",ncoll,")");
my(M=Mat(stats)); for(i=1,#M[,1], print("  ",M[i,1],"  [tested, law holds for all sign-normalized generators, for some]: ",M[i,2]));
}
read("cert.gp");
