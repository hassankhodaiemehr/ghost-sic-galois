default(realprecision,3050);
read("e9h.gp");
{
my(M=36);
foreach([10], p, my(e=e9[p]);
for(k=0,M-1, foreach([1,2],s, my(z=exp(2*Pi*I*k/M)*e^s, t=real(z+1/z), P=algdep(t,48));
  if(abs(subst(P,'x,t))>10^-2500 || 48*log(vecmax(abs(Vec(P))))/log(10)>1800, next);
  my(f=factor(P)[,1], best); for(i=1,#f, if(abs(subst(f[i],'x,t))<10^-2000, best=f[i]));
  print("p=",p," s=",s," k=",k," deg=",poldegree(best)," nonreal=",#select(w->abs(imag(w))>1e-30,polroots(best))," h=10^",round(log(vecmax(abs(Vec(best))))/log(10))))));
}
