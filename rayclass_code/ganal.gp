read("an.gp");
latint(A,B)={ my(K=matkerint(concat(A,-B)), n=#A); mathnf(A*K[1..n,]) }
ghostfields(G,d,Dsf,idxs,maxdeg)={
  my(bnf=bnfinit(y^2-Dsf,1), nf=bnf.nf, rts=nf.roots, idp=if(rts[1]>0,1,2), recs=List(), pr=getlocalprec());
  foreach(idxs, p, my(e=real(G[p]), t=e+1/e, P=algdep(t,maxdeg), sq=0);
    if(abs(subst(P,'x,t))>10^(-pr*0.8) || poldegree(P)*log(vecmax(abs(Vec(P))))/log(10)>pr*0.6,
       e=e^2; t=e+1/e; P=algdep(t,maxdeg); sq=1;
       if(abs(subst(P,'x,t))>10^(-pr*0.8) || poldegree(P)*log(vecmax(abs(Vec(P))))/log(10)>pr*0.6, print("p=",p," not identified up to deg ",maxdeg," (also squared)"); next));
    if(sq, print("p=",p,": using eps^2 (eps itself has larger degree)"));
    my(f=factor(P)[,1], Q); for(i=1,#f, if(abs(subst(f[i],'x,t))<10^(-pr*0.7), Q=f[i]));
    \\ sign-free: e^2
    my(t2=e^2+e^-2, fa=factor(polresultant(Q,'y-('x^2-2),'x)), Q2);
    for(i=1,#fa[,1], my(q=subst(fa[i,1],'y,'x)); if(abs(subst(q,'x,t2))<10^(-pr*0.7), Q2=q));
    my(Pe2=simplify('x^poldegree(Q2)*subst(Q2,'x,'x+1/'x)), Pe=simplify('x^poldegree(Q)*subst(Q,'x,'x+1/'x)));
    my(nr=#select(w->abs(imag(w))>1e-30,polroots(Q)));
    print("p=",p," deg t=",poldegree(Q)," nonreal conj of t: ",nr,"  P_t(0) ",if(polcoef(Q,0)%2,"odd","even"),"  P_t(2)P_t(-2)=",factor(subst(Q,'x,2)*subst(Q,'x,-2)));
    print("   e-field: ",relcond(Pe,e,Dsf));
    print("   e^2-field: ",relcond(Pe2,e^2,Dsf));
    my(fa2=nffactor(nf,Pe2), R=0, best=+oo); for(i=1,#fa2[,1], my(v=abs(subst(subst(lift(fa2[i,1]),'y,rts[idp]),'x,e^2))); if(v<best,best=v;R=lift(fa2[i,1])));
    listput(recs,R));
  if(#recs==0, return);
  my(m=1); for(i=1,#recs, m=idealmul(nf,m,rnfconductor(bnf,recs[i])[1][1]));
  my(bnr=bnrinit(bnf,[m,[1,1]],1), H=matdiagonal(bnr.cyc)); for(i=1,#recs, H=latint(H,rnfnormgroup(bnr,recs[i])));
  my(fm=bnrconductor(bnr,H), sn=select(n->n>1,matsnf(H)));
  print("COMPOSITUM of e^2-fields: [M0:F]=",bnr.no/matdet(H),"  conductor ",idealfactor(nf,fm[1])," inf ",fm[2]," (identity idx ",idp,")  Gal(M0/F)=",sn);
  my(dp=if(d%2,d,2*d)); print("   predicted Cl(d' inf-) = ",bnrinit(bnf,[dp,vector(2,j,j!=idp)]).cyc);
}
