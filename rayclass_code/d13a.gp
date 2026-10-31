default(realprecision,4150);
read("g13.gp");
{
foreach([1,14],p, my(e=real(g13[p]));
 foreach([1,2],s, my(u=e^s, t=u+1/u, P, t0=getabstime());
  P=algdep(t,96);
  my(ok=abs(subst(P,'x,t))<10^-3300 && poldegree(P)*log(vecmax(abs(Vec(P))))/log(10)<2500);
  if(ok, my(f=factor(P)[,1],Q); for(i=1,#f, if(abs(subst(f[i],'x,t))<10^-3000, Q=f[i]));
     print("p=",p," s=",s," deg t=",poldegree(Q)," h=10^",round(log(vecmax(abs(Vec(Q))))/log(10))," time ",(getabstime()-t0)\1000,"s"); write("d13polys.gp","[",p,",",s,",",Q,"]"),
     print("p=",p," s=",s," not found up to 96 (",(getabstime()-t0)\1000,"s)"))));
}
