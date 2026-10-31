read("hs2.gp"); read("gen.gp");
analyzeGhost(E,d)={
  my(N=d^2*(d-3),D=N*(N+4),polys=List(),L=List());
  for(i=1,#E, my(e=E[i][3],t,P,known=0,n,rts,big,nr,q);
    if(abs(abs(e)-1)<1e-50, next);
    listput(L,[0,N,0,i,e]);
    t=e+1/e;
    for(j=1,#polys, if(abs(subst(polys[j],'x,t))<10^(-getlocalprec()*0.7)*vecmax(abs(Vec(polys[j]))), known=1; break));
    if(known,next);
    P=findpol(t,40); if(P==0, print("p=",E[i][1..2]," FAIL"); next); listput(polys,P);
    n=poldegree(P); rts=polroots(P); nr=#select(z->abs(imag(z))>1e-30,rts);
    big=#select(z->abs(imag(z))<=1e-30&&abs(real(z))>2,rts);
    q=signlaw(P,D,t);
    print("p=",E[i][1..2]," deg t=",n," |t|>2:",big," |t|<2:",n-nr-big," nonreal:",nr," sqrt(D) in Q(t), agree, disagree:",q,
      " m=",n/2," P(0)=",polcoef(P,0)," P(2)=",factor(subst(P,'x,2))," P(-2)=",factor(subst(P,'x,-2)), " P:",P));
  Vec(L)
}
