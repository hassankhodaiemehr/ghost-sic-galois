default(realprecision,1150);
read("e9.gp");
{
my(d=9, M=4*d, e=e9[1]);
for(k=0,M-1, my(z=exp(2*Pi*I*k/M)*e, t=real(z+1/z), P=algdep(t,24), ok, f, best);
  if(abs(subst(P,'x,t))>10^-900, next);
  f=factor(P)[,1]; for(i=1,#f, if(abs(subst(f[i],'x,t))<10^-800, best=f[i]));
  my(r=polroots(best)); print("k=",k," deg=",poldegree(best)," nonreal=",#select(w->abs(imag(w))>1e-30,r)," maxcoef~10^",round(log(vecmax(abs(Vec(best))))/log(10))));
}
