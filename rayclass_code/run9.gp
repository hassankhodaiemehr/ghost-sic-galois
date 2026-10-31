default(realprecision,1150);
read("an.gp"); read("e9.gp");
{
for(p=1,4, my(e=e9[p]); for(k=1,2, my(t=real(e^k+e^-k), P=fp(t,80));
  if(P==0, print("p=",p," k=",k," FAIL"),
  print("p=",p," k=",k," deg=",poldegree(P)," in(-2,2): ",#select(z->abs(imag(z))<1e-30 && abs(real(z))<2, polroots(P))," cond: ",relcond(P,t,15)))));
}
