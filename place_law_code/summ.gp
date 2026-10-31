summ(L)={
  my(seen=Map());
  for(i=1,#L, my(r=L[i],P,N=r[2],q,n,a2,am2,rts,big,nr);
    if(abs(abs(r[5])-2)<1e-50, next); P=findpol(r[5],80); if(P==0, print(r[1],"|",N,"|",r[4],"|FAIL"); next);
    if(mapisdefined(seen,Str(P)),next); mapput(seen,Str(P),1);
    n=poldegree(P); rts=polroots(P); nr=#select(z->abs(imag(z))>1e-30,rts);
    big=#select(z->abs(imag(z))<=1e-30&&abs(real(z))>2,rts);
    q=signlaw(P,N*(N+4),r[5]);
    a2=subst(P,'x,2); am2=subst(P,'x,-2);
    print(r[1],"|",N,"|",r[4],"|",n,"|",big,"|",n-nr-big,"|",nr,"|",q,"|",factor(a2),"|",factor(am2)));
}
