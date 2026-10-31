default(realprecision,10200);
read("grefine.gp"); read("gin13.gp");
g=ghostrefine(13,gre13,gim13,10200);
write("g13h.gp","g13=",g,";");
print("refined");
{
foreach([14,1],p, my(e=real(g[p]), u=e^2, t=u+1/u, P, t0=getabstime());
  P=algdep(t,96);
  my(ok=abs(subst(P,'x,t))<10^-8000 && poldegree(P)*log(vecmax(abs(Vec(P))))/log(10)<6500);
  if(ok, my(f=factor(P)[,1],Q); for(i=1,#f, if(abs(subst(f[i],'x,t))<10^-7500, Q=f[i]));
     print("p=",p," s=2 deg t=",poldegree(Q)," h=10^",round(log(vecmax(abs(Vec(Q))))/log(10))," time ",(getabstime()-t0)\1000,"s"); write("d13polys.gp","[",p,",2,",Q,"]"),
     print("p=",p," s=2 not found up to 96 (",(getabstime()-t0)\1000,"s)")));
}
