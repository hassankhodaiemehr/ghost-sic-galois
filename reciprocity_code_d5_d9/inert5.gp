default(realprecision,400);
d=5; D0=3;
L=readvec("g5.gpv");
P=vector(#L,i,[L[i][1],L[i][2]]); V=vector(#L,i,L[i][3]^2);
idx(p)={my(a=p[1]%d,b=p[2]%d); for(i=1,#P, if(P[i][1]==a&&P[i][2]==b, return(i))); 0};
e0=V[idx([0,1])];
{ Pm=algdep(e0,32); fa=factor(Pm)[,1]; for(i=1,#fa, if(abs(subst(fa[i],'x,e0))<1e-300, Pm=fa[i])); }
n=poldegree(Pm);
expr(v)={my(r=lindep(concat([v],vector(n,k,e0^(k-1))))); -Polrev(r[2..n+1])/r[1]};
Ex=vector(#V,i,expr(V[i]));
{
tot=0; good=0; good2=0;
forprime(q=7,3000, if(kronecker(D0,q)!=-1, next); if(poldisc(Pm)%q==0||q==d, next);
  my(g=lift(factormod(Pm,q)[1,1]), red=vector(#P,j,Mod(Mod(1,q)*Ex[j],Mod(1,q)*g)), ok=1, ok2=1);
  for(i=1,#P, my(v=red[i]^(q^2)); if(red[idx(q*P[i])]!=v, ok=0); if(red[idx(-q*P[i])]!=v, ok2=0));
  tot++; good+=ok; good2+=ok2);
print("inert primes: ",tot,"  law G=qI holds: ",good,"   (G=-qI: ",good2,")");
}
