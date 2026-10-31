read("hs2.gp");
condof(x0,D0)={
  my(Px=findpol(x0,64),bnf=bnfinit(y^2-D0,1),nf=bnf.nf,rts=nf.roots,idp=if(rts[1]>0,1,2),fa,R,best=+oo,c);
  fa=nffactor(nf,Px);
  for(i=1,#fa[,1], my(v=abs(subst(subst(lift(fa[i,1]),'y,rts[idp]),'x,x0))); if(v<best,best=v;R=lift(fa[i,1])));
  c=rnfconductor(bnf,R);
  print("  deg=",poldegree(Px)," [M:F]=",poldegree(R,'x)," conductor ",idealfactor(nf,c[1][1])," inf ",c[1][2]," Gal=",c[2].cyc," H=",c[3]);
}
