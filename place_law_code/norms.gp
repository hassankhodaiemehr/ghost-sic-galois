read("hs2.gp");
norms(L,maxd)={
  my(polys=List());
  for(i=1,#L, my(r=L[i],P,N=r[2],known=0);
    if(abs(abs(r[5])-2)<1e-50, next);
    for(j=1,#polys, if(abs(subst(polys[j],'x,r[5]))<10^(-getlocalprec()*0.7)*vecmax(abs(Vec(polys[j]))), known=1; break));
    if(known, next);
    P=findpol(r[5],maxd); if(P==0, next); listput(polys,P);
    print(r[1],"|",N,"|",r[4],"|",poldegree(P),"|",subst(P,'x,2),"|",subst(P,'x,-2),"|",Vec(P)));
}
