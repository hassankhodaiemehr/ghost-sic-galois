read("d13polys.gp");
{
my(Q=Q13, Pu=simplify('x^96*subst(Q,'x,'x+1/'x)), bnf=bnfinit(y^2-35,1), nf=bnf.nf, rts=nf.roots, idp=if(rts[1]>0,1,2), oth=3-idp,
   bm=bnrinit(bnf,[13,vector(2,j,j==oth)]), bp=bnrinit(bnf,[13,vector(2,j,j==idp)]), agree_m=0, dis_m=0, agree_p=0, dis_p=0, n=0, nsplit=0);
print("Cl(13 inf-) ",bm.cyc,"  Cl(13 inf+) ",bp.cyc);
forprime(p=3,6000, if(p==5||p==7||p==13||kronecker(35,p)!=1, next);
  if(poldisc(Pu*Mod(1,p))==0, next);
  my(dec=idealprimedec(nf,p), split=(vecmax(factormod(Pu,p,1)[,1])==1));
  my(tm=prod(k=1,#dec, bnrisprincipal(bm,dec[k],0)==0), tp=prod(k=1,#dec, bnrisprincipal(bp,dec[k],0)==0));
  n++; nsplit+=split; if(split==tm, agree_m++, dis_m++); if(split==tp, agree_p++, dis_p++));
print("split primes p tested: ",n,"  of which P_u splits completely mod p: ",nsplit);
print("prediction (ray class field mod 13*inf_minus): agree ",agree_m," disagree ",dis_m);
print("wrong place (13*inf_plus):                    agree ",agree_p," disagree ",dis_p);
}
