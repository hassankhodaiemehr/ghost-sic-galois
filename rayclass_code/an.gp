read("hs2.gp");
fp(t,maxd)={ my(P,pr=getlocalprec()); for(dd=1,maxd, P=algdep(t,dd);
   if(poldegree(P)==dd && abs(subst(P,'x,t))<10^(-pr*0.8) && dd*log(vecmax(abs(Vec(P)))+1)/log(10)<pr*0.6,
     my(f=factor(P)[,1]); for(i=1,#f, if(abs(subst(f[i],'x,t))<10^(-pr*0.7), return(f[i]))))); 0 }
relcond(P,t0,Dsf)={ \\ field Q(t0) with min poly P ; return [deg, F in Q(t)?, conductor info]
  my(bnf=bnfinit(y^2-Dsf,1),nf=bnf.nf,rts=nf.roots,idp=if(rts[1]>0,1,2),fa,R,best=+oo,c);
  fa=nffactor(nf,P);
  if(#fa[,1]==1, return(["sqrtD not in Q(t)", poldegree(P)]));
  for(i=1,#fa[,1], my(v=abs(subst(subst(lift(fa[i,1]),'y,rts[idp]),'x,t0))); if(v<best,best=v;R=lift(fa[i,1])));
  c=rnfconductor(bnf,R);
  [poldegree(P), poldegree(R,'x), idealfactor(nf,c[1][1]), c[1][2], c[2].cyc, c[3], "identity place idx", idp]
}
