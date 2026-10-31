\\ input: vector of [k,N,sol,u,t] ; precision set by caller
findpol(t,maxd)={
  my(P,pr=getlocalprec());
  for(d=1,maxd, P=algdep(t,d);
    if(poldegree(P)==d && polisirreducible(P) && abs(subst(P,'x,t))<10^(-pr*0.8) && d*log(vecmax(abs(Vec(P)))+1)/log(10) < pr*0.55, return(P)));
  0
}
\\ for field Q(t) given by P, test sqrt(q) in Q(t) and sign law
signlaw(P,q,t0)={
  my(C=polcompositum(P,'x^2-q,1),c,R,rts,a,b,agree=0,dis=0,n,o,best=10^9,bi);
  if(#C!=2, return([0,0,0]));
  c=C[1]; R=c[1]; a=c[2]; b=c[3]; rts=polroots(R); n=#rts;
  for(i=1,n, my(dd=abs(subst(lift(a),'x,rts[i])-t0)); if(dd<best,best=dd;bi=i));
  o=sign(real(subst(lift(b),'x,rts[bi])));
  for(i=1,n, my(tv=real(subst(lift(a),'x,rts[i])), sv=real(subst(lift(b),'x,rts[i])));
     if(sign(tv^2-4)==o*sign(sv), agree++, dis++));
  [1,agree,dis]
}
analyze(L)={
  my(seen=Map(),res=List());
  for(i=1,#L, my(r=L[i],P,key,rts,n,big,small,nonreal,q1,q2,q3,N=r[2]);
    if(abs(abs(r[5])-2)<1e-50, next); P=findpol(r[5],80);
    if(P==0, print("k=",r[1]," N=",N," sol=",r[3]," u=",r[4],": no poly"); next);
    key=Str(P); if(mapisdefined(seen,key), next); mapput(seen,key,1);
    rts=polroots(P); n=poldegree(P);
    nonreal=#select(z->abs(imag(z))>1e-30,rts);
    big=#select(z->abs(imag(z))<=1e-30 && abs(real(z))>2,rts); small=n-nonreal-big;
    q1=signlaw(P,N+4,r[5]); q2=signlaw(P,N,r[5]); q3=signlaw(P,N*(N+4),r[5]);
    print("k=",r[1]," N=",N," sol=",r[3]," u=",r[4]," deg=",n," |t|>2:",big," |t|<2:",small," nonreal:",nonreal,
          " sqrt(N+4):",q1," sqrt(N):",q2," sqrt(N(N+4)):",q3," P2P-2:",sign(subst(P,'x,2)*subst(P,'x,-2))));
}
