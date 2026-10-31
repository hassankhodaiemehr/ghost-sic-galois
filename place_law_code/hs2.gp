findpol(t,maxd)={
  my(P,pr=getlocalprec());
  forstep(d=2,maxd,2, P=algdep(t,d);
    if(poldegree(P)==d && abs(subst(P,'x,t))<10^(-pr*0.8) && d*log(vecmax(abs(Vec(P)))+1)/log(10) < pr*0.55,
       my(f=factor(P)[,1]); for(i=1,#f, if(abs(subst(f[i],'x,t))<10^(-pr*0.7), return(f[i])))));
  0
}
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
summ(L,maxd)={
  my(polys=List());
  for(i=1,#L, my(r=L[i],P=0,N=r[2],q,n,rts,big,nr,known=0);
    if(abs(abs(r[5])-2)<1e-50, next);
    for(j=1,#polys, if(abs(subst(polys[j],'x,r[5]))<10^(-getlocalprec()*0.7)*vecmax(abs(Vec(polys[j]))), known=1; break));
    if(known, next);
    P=findpol(r[5],maxd); if(P==0, print(r[1],"|",N,"|",r[4],"|FAIL>",maxd); next);
    listput(polys,P);
    n=poldegree(P); rts=polroots(P); nr=#select(z->abs(imag(z))>1e-30,rts);
    big=#select(z->abs(imag(z))<=1e-30&&abs(real(z))>2,rts);
    q=signlaw(P,N*(N+4),r[5]);
    print(r[1],"|",N,"|",r[4],"|",n,"|",big,"|",n-nr-big,"|",nr,"|",q[1],"|",q[2],"|",q[3],"|",sign(subst(P,'x,2)*subst(P,'x,-2))));
}
