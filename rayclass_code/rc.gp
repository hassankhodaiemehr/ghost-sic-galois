{
for(i=1,5, my(d=[4,5,9,13,17][i], Dsf=core((d-3)*(d+1)), dp=if(d%2,d,2*d), bnf, r, idm, bnr, bnrb, bnr0);
  bnf=bnfinit(y^2-Dsf,1); r=bnf.nf.roots; idm=if(r[1]<0,1,2);
  bnr=bnrinit(bnf,[dp,vector(2,j,j==idm)]); bnrb=bnrinit(bnf,[dp,[1,1]]); bnr0=bnrinit(bnf,dp);
  print("d=",d," F=Q(sqrt ",Dsf,") h=",bnf.no," d'=",dp,"  |Cl(d' inf-)|=",bnr.no," cyc ",bnr.cyc," |Cl(d' inf+inf-)|=",bnrb.no," |Cl(d')|=",bnr0.no));
}
