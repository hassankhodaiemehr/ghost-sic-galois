default(parisize,"1G");
{
my(W=250, tot=vector(8), hit=vector(8), L=List());
for(d=4,2000, my(Dsf=core((d-3)*(d+1)), dp=if(d%2,d,2*d), bnf, r, idm, b, b0, c, rk2, w=(d-1)\W+1);
  if(Dsf==1, next);
  bnf=bnfinit(y^2-Dsf,1); r=bnf.nf.roots; idm=if(r[1]<0,1,2);
  b=bnrinit(bnf,[dp,vector(2,j,j==idm)]); b0=bnrinit(bnf,dp); c=b.cyc;
  rk2=#select(n->n%2==0,c);
  tot[w]++; if(rk2==1 && b.no==2*b0.no, hit[w]++; listput(L,[d,omega(Dsf),isprime(d)])));
for(i=1,8, print("d in [",(i-1)*W+1,",",i*W,"]: ",hit[i],"/",tot[i]));
print(L);
}
