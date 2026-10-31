\\ L = [[k,N,sol,u,x],...]
\\ For each u0: (gen) x_v in Q(x_u0) for all v ; (rec) P_{x_u0} reciprocal ;
\\ (R2') kappa: x_u0 -> 1/x_u0 maps every x_v to 1/x_v (exact check mod P_{x_u0}).
expr(P0,x0,Pv,xv)={ \\ return f in Q[x] with f(x0)=xv, or 0 if xv not in Q(x0)
  my(C=polcompositum(P0,Pv,1),best=+oo,bc=0);
  for(c=1,#C, my(R=C[c][1],k=C[c][4],th=xv+k*x0,v);
    v=log(abs(subst(R,'x,th))+10^(-getlocalprec()))-log(vecmax(abs(Vec(R)))*(1+abs(th))^poldegree(R));
    if(v<best,best=v;bc=c));
  my(R=C[bc][1],a=C[bc][2],b=C[bc][3]);
  if(poldegree(R)!=poldegree(P0), return(0));
  my(th=modreverse(a)); \\ theta as polmod in x0 (mod charpoly of a = P0)
  lift(subst(lift(b),'x,th))
}
genanalyze(L)={
  my(n=#L,P=vector(n));
  for(i=1,n, P[i]=if(abs(abs(L[i][5])-1)<1e-50,0,findpol(L[i][5],120)));
  for(i=1,n, if(P[i]==0,next);
    my(P0=P[i],gen=1,r2=1,y=Mod('x,P0),rec=(polrecip(P0)==P0||polrecip(P0)==-P0));
    for(j=1,n, if(j==i||P[j]==0,next);
      my(f=expr(P0,L[i][5],P[j],L[j][5]));
      if(f==0, gen=0; r2=-1; break);
      if(abs(subst(f,'x,L[i][5])-L[j][5])>1e-100, print("  WARNING mismatch v=",L[j][4]));
      if(subst(f,'x,y)*subst(f,'x,y^-1)!=1, r2=0));
    print("k=",L[i][1]," N=",L[i][2]," u0=",L[i][4]," deg x=",poldegree(P0)," reciprocal=",rec," generates=",gen," kappa inverts all x_v=",r2,
          " P_t(2)P_t(-2) sign=", if(rec, my(Pt=findpol(L[i][5]+1/L[i][5],80)); sign(subst(Pt,'x,2)*subst(Pt,'x,-2)), "na")));
}
