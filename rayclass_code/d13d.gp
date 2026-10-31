default(realprecision,20500);
read("an.gp"); read("d13polys.gp"); read("g13x.gp");
{
my(e=real(g13[14]), u=e^2, t=u+1/u, Q=Q13, r=polroots(Q));
print("deg ",poldegree(Q)," nonreal ",#select(w->abs(imag(w))>1e-30,r)," in(-2,2) ",#select(w->abs(imag(w))<1e-30&&abs(real(w))<2,r)," P(0) ",if(polcoef(Q,0)%2,"odd","even"));
my(t0=getabstime());
print("t-field: ", relcond(Q,t,35), "  (",(getabstime()-t0)\1000,"s)");
my(Pu=simplify('x^96*subst(Q,'x,'x+1/'x)));
t0=getabstime();
print("u=eps^2 field: ", relcond(Pu,u,35), "  (",(getabstime()-t0)\1000,"s)");
my(bnf=bnfinit(y^2-35,1)); print("Cl(13 inf-) predicted: ",bnrinit(bnf,[13,[1,0]]).cyc," / [0,1]: ",bnrinit(bnf,[13,[0,1]]).cyc, " Cl(13): ",bnrinit(bnf,13).cyc);
my(a=subst(Q,'x,2), b=subst(Q,'x,-2)); print("P(2)P(-2) digits ",#Str(abs(a*b)));
print("factor P(2): ",factor(a,10^7)); print("factor P(-2): ",factor(b,10^7));
}
