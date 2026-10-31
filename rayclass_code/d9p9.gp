default(realprecision,1150);
read("e9.gp");
{
my(d=9, M=36);
foreach([10,18], p, my(e=e9[p]);
for(k=0,M-1, my(z=exp(2*Pi*I*k/M)*e^2, t=real(z+1/z), P=algdep(t,36));
  if(abs(subst(P,'x,t))>10^-900 || 36*log(vecmax(abs(Vec(P))))/log(10)>600, next);
  my(f=factor(P)[,1], best); for(i=1,#f, if(abs(subst(f[i],'x,t))<10^-800, best=f[i]));
  print("p=",p," k=",k,"/36 deg=",poldegree(best)," nonreal=",#select(w->abs(imag(w))>1e-30,polroots(best)))));
}
