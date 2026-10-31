read("../rayclass_code/g9.gp"); read("aut9.gp");
default(realprecision,2000); TOL=1e-1500;
d=9; D0=15;
e=vector(80,i,real(g9[i])^2);
P=List(); for(a=0,8, for(b=0,8, if(a||b, listput(P,[a,b])))); P=Vec(P);
idx(p)={my(a=p[1]%d,b=p[2]%d); for(i=1,#P, if(P[i][1]==a&&P[i][2]==b, return(i))); 0};
mp(v,n)={my(Q=algdep(v,n),f=factor(Q)[,1]); for(k=1,#f, if(abs(subst(f[k],'x,v))<TOL, return(f[k]))); 0};
e0=e[idx([1,0])]; Pm=Pm9; Pmt=subst(Pm,'x,'t);
read("ex9_all.gp"); Ex=Ex9;
\\ numeric action of each automorphism on the overlaps of this ghost: perm or 0 if leaves the ghost
{ perms=vector(#Aut9); for(k=1,#Aut9, my(t=subst(Aut9[k],'x,e0), pm=vector(80)); 
   for(i=1,80, my(v=subst(Ex[i],'x,t)); for(j=1,80, if(abs(v-e[j])<TOL, pm[i]=j; break)));
   perms[k]=pm);
  print("automorphisms preserving this ghost: ", #select(p->vecmin(p)>0, perms)," of ",#Aut9);
  print("automorphisms moving all overlaps outside: ", #select(p->vecmax(p)==0, perms)); }
{
S=List();
forvec(m=vector(4,i,[0,d-1]), my(G=[m[1],m[2];m[3],m[4]],ok=1);
  if(gcd(matdet(G),d)!=1,next);
  for(i=1,#P, my(q=lift(Mod(G,d)*P[i]~)); if(abs(e[idx(q)]-e[i])>TOL, ok=0;break));
  if(ok, listput(S,G)));
Lz=0; for(i=1,#S, my(G=S[i]); if(trace(G)%d==(d-1)%d && matdet(G)%d==1 && G!=matid(2), Lz=G)); print("L = ",Lz);
}
Xs=[Lz-4*matid(2), lift(Mod(Lz,d)^-1)-4*matid(2)];
print("X^2 = ",lift(Mod(Xs[1],d)^2)," ",lift(Mod(Xs[2],d)^2));
bnf=bnfinit(y^2-D0,1); dens=lcm(apply(c->denominator(content(c)),concat(Ex,Vec(Aut9)))); discP=poldisc(Pm);
\\ sqrt(15) as element: via automorphism-free nfroots
{ my(R=nfroots(Pmt,'x^2-D0)); sq=0; for(k=1,#R, my(r=subst(lift(R[k]),'t,'x)); if(abs(subst(r,'x,e0)-sqrt(D0))<TOL, sq=r)); print("sq ok ",sq!=0); }
{
nofrob=0; stats=Map(); total=0; nonpr=0; nonprPres=0; prPres=0; inert=Map(); ni=0; niPres=0;
forprime(q=7,4000, if(D0%q==0||dens%q==0||discP%q==0, next);
  my(g=lift(factormod(Pm,q)[1,1]), md=Mod(1,q)*g, fr=0);
  my(tq=Mod(Mod(1,q)*'x,md)^if(kronecker(D0,q)==-1,q^2,q));
  for(k=1,#Aut9, if(Mod(Mod(1,q)*Aut9[k],md)==tq, fr=k; break));
  if(!fr, nofrob++; next);
  my(pm=perms[fr], pres=(vecmin(pm)>0));
  if(kronecker(D0,q)==-1, ni++; niPres+=pres; if(pres,
      foreach([1,-1],s, my(G=Mod(s*q*matid(2),d), ok=1); for(i=1,80, if(idx(lift(G*P[i]~))!=pm[i] && abs(e[idx(lift(G*P[i]~))]-e[pm[i]])>TOL, ok=0; break)); mapput(inert,s,if(mapisdefined(inert,s),mapget(inert,s),0)+ok)));
      next);
  my(rv=lift(lift(Mod(Mod(1,q)*sq,md)))); if(poldegree(rv)>0, print("sqrt not in Fq at ",q); next); rv=polcoef(rv,0);
  my(pr=idealprimedec(bnf,q), pid=0); for(k=1,#pr, if(nfeltval(bnf,'y-rv,pr[k])>0, pid=pr[k]));
  my(ip=bnfisprincipal(bnf,pid));
  if(ip[1]!=0, nonpr++; nonprPres+=pres; next);
  prPres+=pres; if(!pres, next);
  my(pi=nfbasistoalg(bnf,ip[2]), a=polcoef(lift(pi),0), b=polcoef(lift(pi),1));
  if(denominator(a)>1||denominator(b)>1, print("nonintegral gen"); next);
  total++;
  foreach([1,2],xi, foreach([1,-1],sgn, foreach([1,-1],inv, foreach(["+","-"],place,
     my(A=a,B=b, s=if(place=="+", sign(A+B*sqrt(D0)), sign(A-B*sqrt(D0)))); A*=s; B*=s;
     my(G=Mod(A*matid(2)+B*Xs[xi],d)); if(inv==-1, G=G^-1); if(sgn==-1, G=G~);
     my(ok=1); for(i=1,80, my(j=idx(lift(G*P[i]~))); if(abs(e[j]-e[pm[i]])>TOL, ok=0; break));
     my(key=Str(xi,sgn,inv,place)); mapput(stats,key, if(mapisdefined(stats,key),mapget(stats,key),0)+ok))))));
print("split principal (ghost-preserving) tested: ",total,"; principal primes total preserve: ",prPres);
print("nonprincipal split primes: ",nonpr,", of which ghost-preserving: ",nonprPres);
print("no frob found: ",nofrob); print("inert: ",ni," ghost-preserving ",niPres," law qI/-qI: ",Mat(inert));
my(M=Mat(stats)); for(i=1,#M[,1], print(M[i,1]," -> ",M[i,2]));
}
