{
my(cnt=0, cyc2=0, list=List());
for(d=4,150, my(Dsf=core((d-3)*(d+1)), dp=if(d%2,d,2*d), bnf, r, idm, b, b0, c, v2, rk2);
  if(Dsf==1, next);
  bnf=bnfinit(y^2-Dsf,1); r=bnf.nf.roots; idm=if(r[1]<0,1,2);
  b=bnrinit(bnf,[dp,vector(2,j,j==idm)]); b0=bnrinit(bnf,dp); c=b.cyc;
  rk2=#select(n->n%2==0,c);
  cnt++; if(rk2==1 && b.no==2*b0.no, cyc2++; listput(list,d));
  if(d<=40 || rk2==1, print("d=",d," F=Q(sqrt ",Dsf,") Cl(d' inf-)=",c,"  |.|/|Cl(d')|=",b.no/b0.no,"  2-rank=",rk2)));
print("dims 4..150 tested: ",cnt,"; 2-Sylow cyclic and inf- in conductor: ",#list);
print(Vec(list));
}
