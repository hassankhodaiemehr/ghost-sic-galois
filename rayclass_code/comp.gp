read("an.gp");
\\ for list of SIC overlaps E (unimodular), dimension d: find twist, sign-free unit field per overlap, compositum data
latint(A,B)={ my(K=matkerint(concat(A,-B)), n=#A); mathnf(A*K[1..n,]) }
unitdata(E,d,Dsf,maxdeg)={
  my(M=4*d, bnf=bnfinit(y^2-Dsf,1), nf=bnf.nf, rts=nf.roots, idp=if(rts[1]>0,1,2), seen=List(), recs=List());
  for(p=1,#E, my(e=E[p], got=0);
    for(k=0,M-1, if(got,break); my(z=exp(2*Pi*I*k/M)*e, t=real(z+1/z), P=algdep(t,maxdeg));
      if(abs(subst(P,'x,t))>10^(-getlocalprec()*0.8), next);
      my(f=factor(P)[,1], Q=0); for(i=1,#f, if(abs(subst(f[i],'x,t))<10^(-getlocalprec()*0.7), Q=f[i]));
      if(Q==0, next);
      if(#select(w->abs(imag(w))>1e-30, polroots(Q))>0, next);   \\ require totally real t
      got=1;
      \\ sign-free: t2=u^2+u^-2
      my(t2=t^2-2, fa=factor(polresultant(Q,'y-('x^2-2),'x)), Q2=0);
      for(i=1,#fa[,1], my(q=subst(fa[i,1],'y,'x)); if(abs(subst(q,'x,t2))<10^(-getlocalprec()*0.7), Q2=q));
      my(Pu2=simplify('x^poldegree(Q2)*subst(Q2,'x,'x+1/'x)));
      my(key=Str(Pu2)); if(setsearch(Set(Vec(seen)),key), next); listput(seen,key);
      \\ relative polynomial over F of Q(u^2)
      my(fa2=nffactor(nf,Pu2), R=0, best=+oo, u2=z^2);
      for(i=1,#fa2[,1], my(v=abs(subst(subst(lift(fa2[i,1]),'y,rts[idp]),'x,u2))); if(v<best,best=v;R=lift(fa2[i,1])));
      listput(recs,[p,k,poldegree(Q),poldegree(Pu2),R]));
    if(!got, print("p=",p," no totally real twist found up to deg ",maxdeg)));
  print("distinct sign-free unit fields: ",#recs);
  \\ common modulus: lcm of conductors
  my(conds=List(), m=1, inf=[0,0]);
  for(i=1,#recs, my(c=rnfconductor(bnf,recs[i][5])); listput(conds,c); m=idealmul(nf,m,c[1][1]); inf=[max(inf[1],c[1][2][1]),max(inf[2],c[1][2][2])];
     print("  overlap p=",recs[i][1]," twist k=",recs[i][2]," deg t=",recs[i][3]," [Q(u^2):Q]=",recs[i][4]," conductor ",idealfactor(nf,c[1][1])," inf ",c[1][2]," Gal ",c[2].cyc));
  my(bnr=bnrinit(bnf,[m,[1,1]],1), H=matdiagonal(bnr.cyc));
  for(i=1,#recs, H=latint(H, rnfnormgroup(bnr,recs[i][5])));
  my(G=bnrinit(bnf,[m,[1,1]],1)); 
  print("compositum: [M0:F] = ", bnr.no/matdet(H), "  (Cl_m order ",bnr.no,")");
  my(fm=bnrconductor(bnr,H,1)); print("compositum conductor: ",idealfactor(nf,fm[1][1])," inf ",fm[1][2], "  Gal(M0/F) = ", fm[2].cyc);
  my(dp=if(d%2,d,2*d)); for(s=0,1, my(b=bnrinit(bnf,[dp,if(s,[1,1],vector(2,j,j!=idp))])); print("   ray class group mod d'",if(s," inf+inf-"," inf(non-identity)"),": ",b.cyc));
}
