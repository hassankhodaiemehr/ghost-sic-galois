default(realprecision,400);
d=5; D0=3;
L=readvec("g5.gpv");
P=vector(#L,i,[L[i][1],L[i][2]]); V=vector(#L,i,L[i][3]^2);
idx(p)={my(a=p[1]%d,b=p[2]%d); for(i=1,#P, if(P[i][1]==a&&P[i][2]==b, return(i))); 0};
{
S=List();
forvec(m=vector(4,i,[0,d-1]), my(G=[m[1],m[2];m[3],m[4]],ok=1);
  if(matdet(Mod(G,d))==0,next);
  for(i=1,#P, my(q=lift(Mod(G,d)*P[i]~)); if(abs(V[idx(q)]-V[i])>1e-300, ok=0;break));
  if(ok, listput(S,G)));
print("stabilizer S size ",#S); for(i=1,#S, print("  ",S[i]));
}
e0=V[idx([0,1])];
{ Pm=algdep(e0,32); fa=factor(Pm)[,1]; for(i=1,#fa, if(abs(subst(fa[i],'x,e0))<1e-300, Pm=fa[i])); }
n=poldegree(Pm); print("minpoly deg ",n, "  ", Pm);
expr(v)={my(r=lindep(concat([v],vector(n,k,e0^(k-1))))); -Polrev(r[2..n+1])/r[1]};
Ex=vector(#V,i,expr(V[i]));
sq=expr(sqrt(D0));
{ for(i=1,#V, if(abs(subst(Ex[i],'x,e0)-V[i])>1e-300, print("bad expr ",i))); }
d=5; D0=3; Lz=[1,1;2,3];
idx(p)={my(a=p[1]%d,b=p[2]%d); for(i=1,#P, if(P[i][1]==a&&P[i][2]==b, return(i))); 0};
Xs=[Lz-2*matid(2), Lz^2-2*matid(2)];
\\ value classes: index -> class id (equal values)
{
stats=Map(); total=0; lin=0;
forprime(q=7,3000, if(kronecker(D0,q)!=1, next); if(poldisc(Pm)%q==0, next);
  my(fa=factormod(Pm,q)[,1]);
  for(t=1,#fa, my(g=lift(fa[t]), f=poldegree(g), r, perm=vector(#P));
    if(f==1, next); \\ Frobenius trivial on this prime
    r=lift(Mod(subst(sq,'x,Mod(Mod(x,q),Mod(g,q))),1)); 
    my(rv=lift(lift(Mod(subst(lift(sq),'x,Mod('x,Mod(g,q)) ),1))));
    rv=subst(lift(sq),'x,Mod(Mod(1,q)*'x,Mod(1,q)*g)); rv=lift(lift(rv)); if(poldegree(rv)>0, print("sqrt3 not in Fq?"); next); rv=polcoef(rv,0);
    for(i=1,#P, my(vq=Mod(Mod(1,q)*subst(Ex[i],'x,'x),Mod(1,q)*g)^q, found=0);
      for(j=1,#P, if(Mod(Mod(1,q)*Ex[j],Mod(1,q)*g)==vq, perm[i]=j; found=1; break)); if(!found, perm[i]=0));
    if(vecmin(perm)==0, print("q=",q," perm incomplete"); next);
    total++;
    \\ find generator pi=a+b*sqrt3 of prime (q, sqrt3 - rv): need a + b*rv = 0 mod q, a^2-3b^2 = +-q
    my(pi=0); for(b=1,400, my(a0=lift(Mod(-b*rv,q)), A); foreach([a0,a0-q],A, if(abs(A^2-3*b^2)==q, pi=[A,b]; break)); if(pi,break));
    if(!pi, print("no generator found q=",q); next);
    \\ test all conventions
    foreach([1,2],xi, foreach([1,-1],sgn, foreach([1,-1],inv, foreach(["+","-"],place,
       my(a=pi[1],b=pi[2]);
       my(vplus=a+b*sqrt(3.), vminus=a-b*sqrt(3.));
       my(s=if(place=="+", sign(vplus), sign(vminus))); a*=s; b*=s;  \\ make pi positive at chosen place
       my(G=Mod(a*matid(2)+b*Xs[xi],d)); if(inv==-1, G=G^-1); if(sgn==-1, G=G~);
       my(ok=1); for(i=1,#P, my(q2=lift(G*P[i]~)); if(abs(1)&& !(subst(Ex[idx(q2)],'x,'x)==Ex[perm[i]] || idx(q2)==perm[i]), ok=0; break));
       my(key=Str(xi,sgn,inv,place)); mapput(stats,key, if(mapisdefined(stats,key),mapget(stats,key),0)+ok)))))));
print("primes tested (nontrivial Frobenius): ",total);
my(M=Mat(stats)); for(i=1,#M[,1], print(M[i,1]," -> ",M[i,2]));
}
