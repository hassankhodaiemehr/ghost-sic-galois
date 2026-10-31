default(realprecision,1150);
read("an.gp"); read("e9.gp");
{
my(e=e9[1], t=real(e^2+e^-2), P=fp(t,40), r=polroots(P));
print("t deg ",poldegree(P)," nonreal: ",#select(z->abs(imag(z))>1e-30,r)," |t|<2 real: ",#select(z->abs(imag(z))<1e-30&&abs(real(z))<2,r)," |t|>2: ",#select(z->abs(imag(z))<1e-30&&abs(real(z))>2,r));
print("P(2)=",factor(subst(P,'x,2))," P(-2)=",factor(subst(P,'x,-2)));
my(u=e^2, Q=fp(u,60));
print("u=e^2 deg ",poldegree(Q), " cond ", relcond(Q,u,15));
}
