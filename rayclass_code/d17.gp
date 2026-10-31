default(realprecision,20500);
read("grefine.gp"); read("gin17.gp");
t0=getabstime(); g=ghostrefine(17,gre17,gim17,20500); print("refined in ",(getabstime()-t0)\1000,"s; max imag ",vecmax(apply(z->abs(imag(z)),g)));
write("g17x.gp","g17=",g,";");
{
my(found=0);
foreach([18,19,1],p, if(found,break); my(e=real(g[p]), u=e^2, t=u+1/u, P, t0=getabstime());
  P=algdep(t,96);
  my(ok=abs(subst(P,'x,t))<10^-16000 && poldegree(P)*log(vecmax(abs(Vec(P))))/log(10)<13000);
  if(ok, my(f=factor(P)[,1],Q); for(i=1,#f, if(abs(subst(f[i],'x,t))<10^-15000, Q=f[i]));
     print("p=",p," deg t=",poldegree(Q)," h=10^",round(log(vecmax(abs(Vec(Q))))/log(10))," time ",(getabstime()-t0)\1000,"s"); write("d17polys.gp","Q17=",Q,";\np17=",p,";"); found=1,
     print("p=",p," not found up to 96 (",(getabstime()-t0)\1000,"s)")));
}
