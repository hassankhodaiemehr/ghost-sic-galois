\\ Place-law sign check for d=9: every conjugate of t=e+1/e (e = squared ghost overlap)
\\ satisfies sign(t^2-4) = sign of the matching conjugate of sqrt(15).
read("hs2.gp");
default(realprecision,2000);
read("g9.gp");
{
my(e=real(g9[10])^2, t=e+1/e, P=algdep(t,36), f=factor(P)[,1], Q);
for(i=1,#f, if(abs(subst(f[i],'x,t))<10^-1500, Q=f[i]));
print("d=9 deg ",poldegree(Q)," signlaw [ok,agree,disagree]: ",signlaw(Q,15,t));
}
