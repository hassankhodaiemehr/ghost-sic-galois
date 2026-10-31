read("hs2.gp");
\\ x0: high-precision real value of the unit at the identity embedding; D0 squarefree part of D
cft(x0,D0,pmax)={
  my(Px=findpol(x0,64),n=poldegree(Px),bnf,nf,fa,R,best=+oo,auts,cond,bnr,H,rts,idp,fin,inf,beta,cb,stats=Map());
  print("deg x = ",n);
  bnf=bnfinit(y^2-D0,1); nf=bnf.nf;
  \\ identity embedding of F: sqrt(D0)>0 ; find index of real place with y>0
  rts=nf.roots; idp=if(rts[1]>0,1,2);
  print("real places of F (values of y): ",rts,"  identity place index = ",idp);
  fa=nffactor(nf,subst(Px,'x,'x)); 
  \\ choose factor vanishing at x0 under identity embedding
  for(i=1,#fa[,1], my(f=fa[i,1],v); v=abs(subst(subst(lift(f),'y,rts[idp]),'x,x0)); if(v<best,best=v;R=lift(f)));
  print("relative degree [M:F] = ",poldegree(R,'x));
  auts=nfgaloisconj(Px); print("#Aut(M/Q) = ",#auts);
  \\ commutativity of automorphisms
  my(ab=1); for(i=1,#auts, for(j=1,#auts, if(subst(auts[i],'x,Mod(auts[j],Px))!=subst(auts[j],'x,Mod(auts[i],Px)), ab=0)));
  print("Aut(M) abelian: ",ab);
  cond=rnfconductor(bnf,R); fin=cond[1][1]; inf=cond[1][2];
  print("conductor finite part: ",idealfactor(nf,fin), "  norm ",idealnorm(nf,fin),"  infinite part ",inf," (1 = ramified), identity place ",idp);
  bnr=cond[2]; H=cond[3];
  \\ c_-: class of (beta), beta = 1 mod fin, negative at non-identity place, positive at identity place
  my(oth=3-idp, z=nfbasistoalg(nf,idealtwoelt(nf,fin)[2]), g, found=0);
  for(a=-60,60, for(b=-60,60, if(found,break); g=1+ (a+b*'y)*idealnorm(nf,fin);
     my(s1=subst(lift(g),'y,rts[idp]), s2=subst(lift(g),'y,rts[oth]));
     if(s1>0 && s2<0 && norm(Mod(g,y^2-D0))!=0, beta=g; found=1)));
  cb=bnrisprincipal(bnr,idealhnf(nf,beta),0);
  print("beta = ",beta,"  class of c_- in ray group: ",cb,"   bnr.cyc=",bnr.cyc,"  H=",H);
  \\ quotient map: class mod H ; we compare classes via H: two classes equal in Gal(M/F) iff difference in image of H
  my(Hm=matconcat([H,matdiagonal(bnr.cyc)]), same(u,v)=my(w=u-v); type(matsolvemod(H,bnr.cyc~,w))=="t_COL");
  my(bp,f2=0); for(a=-60,60, for(b=-60,60, if(f2,break); my(g=1+(a+b*'y)*idealnorm(nf,fin), s1=subst(lift(g),'y,rts[idp]), s2=subst(lift(g),'y,rts[oth])); if(s1<0 && s2>0, bp=g; f2=1)));
  print("wrong-place element beta_+ (neg at identity place, pos at other): class ",bnrisprincipal(bnr,idealhnf(nf,bp),0), " -> Galois element is ", if(same(bnrisprincipal(bnr,idealhnf(nf,bp),0),vector(#bnr.cyc,i,0)~),"trivial","nontrivial"));
  my(nk=0,nid=0,nbad=0,nother=0, tot=0);
  forprime(p=3,pmax, my(dec=idealprimedec(nf,p));
    for(k=1,#dec, my(P=dec[k]); if(P.f!=1,next); if(idealval(nf,fin,P)>0,next);
      my(modP=nfmodprinit(nf,P), Rp, fr, cl, typ);
      Rp=Pol(apply(c->nfmodpr(nf,c,modP),Vec(R)),'x); if(poldisc(Rp)==0, next);
      fr=factor(Rp); 
      typ = if(vecmax(apply(poldegree,fr[,1]))==1, "split",
             if(vecmax(apply(poldegree,fr[,1]))==2 && vecmin(apply(poldegree,fr[,1]))==2 && prod(i=1,#fr[,1],polcoef(fr[i,1],0))==1 && #select(f->polcoef(f,0)!=1,fr[,1])==0, "kappa","other"));
      cl=bnrisprincipal(bnr,P,0);
      my(pred = if(same(cl,vector(#cl,i,0)~),"split", if(same(cl,cb),"kappa","other")));
      tot++; mapput(stats,Str(pred,"->",typ), if(mapisdefined(stats,Str(pred,"->",typ)),mapget(stats,Str(pred,"->",typ)),0)+1)));
  print("degree-1 primes tested: ",tot);
  print("predicted (from ray class) -> observed Frobenius on x: ",Mat(stats));
}
