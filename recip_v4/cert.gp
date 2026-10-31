\\ appended after recip.gp: exact class-field certificate and ray-class coverage
{
my(rts=bnf.roots, im=if(rts[1]<0,1,2), ip=3-im, arch=vector(2,j,j==im));
print("PARI place ",im," is inf- (sqrt D0 -> ",sign(rts[im]),")");
\\ relative polynomial of e0 over F, chosen at the ghost embedding (y -> +sqrt D0)
my(fr=nffactor(bnf, Pm)[,1], rel=0);
for(i=1,#fr, my(g=fr[i]); if(abs(subst(subst(lift(g),y,sqrt(D0)),'x,e0))<eps, rel=g));
print("relative degree [K:F] = ",poldegree(rel));
my(cond=rnfconductor(bnf, rel)); 
print("conductor of K/F: finite part ",idealfactor(bnf,cond[1][1])," ; archimedean ",cond[1][2]);
bnrM=bnrinit(bnf,[d,arch],1); print("Cl(d inf-) = ",bnrM.cyc);
\\ coverage: ray classes hit by tested split/inert primes, with the law verified (key Lam^-1 pos@inf-)
my(seen=Map());
forprime(q=3,QMAX, if(d%q==0||Dl%q==0||poldisc(Pm)%q==0||den%q==0, next);
  foreach(idealprimedec(bnf,q),pr, my(c=bnrisprincipal(bnrM,pr,0)); if(!mapisdefined(seen,c), mapput(seen,c,q))));
print("ray classes hit by tested primes: ",#Mat(seen)[,1]," of ",vecprod(bnrM.cyc));
}
