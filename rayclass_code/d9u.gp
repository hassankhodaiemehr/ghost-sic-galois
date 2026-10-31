default(realprecision,1150);
read("an.gp"); read("e9.gp");
{
my(d=9, e=e9[1], z=exp(2*Pi*I*14/36)*e, t=real(z+1/z), P, f, Pu, Pt2, Pu2);
P=algdep(t,24); f=factor(P)[,1]; for(i=1,#f, if(abs(subst(f[i],'x,t))<10^-800, P=f[i]));
print("P_t = ",P);
print("t-field: ", relcond(P,t,15));
Pu=numerator(substpol(subst(P,'x,('x^2+1)/'x)*'x^poldegree(P),'x,'x)); 
Pu='x^poldegree(P)*subst(P,'x,'x+1/'x); Pu=simplify(Pu);
print("u poly irreducible: ",polisirreducible(Pu), " deg ",poldegree(Pu));
print("u-field: ", relcond(Pu,z,15));
Pt2=polresultant(P,'y-('x^2-2),'x); Pt2=subst(Pt2,'y,'x); Pt2=factor(Pt2)[1,1];
my(t2=t^2-2, fa=factor(polresultant(P,'y-('x^2-2),'x))); 
for(i=1,#fa[,1], my(q=subst(fa[i,1],'y,'x)); if(abs(subst(q,'x,t2))<10^-800, Pt2=q));
Pu2=simplify('x^poldegree(Pt2)*subst(Pt2,'x,'x+1/'x));
print("u^2-field: deg ",poldegree(Pu2)," ", relcond(Pu2,z^2,15));
}
