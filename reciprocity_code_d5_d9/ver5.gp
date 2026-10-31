default(realprecision,400);
d=5; D0=3;
L=readvec("g5.gpv");
V=vector(#L,i,L[i][3]^2);
e0=V[1];
{ Pm=algdep(e0,32); fa=factor(Pm)[,1]; for(i=1,#fa, if(abs(subst(fa[i],'x,e0))<1e-300, Pm=fa[i])); }
n=poldegree(Pm);
expr(v)={my(r=lindep(concat([v],vector(n,k,e0^(k-1))))); -Polrev(r[2..n+1])/r[1]};
sq=expr(sqrt(D0)); print("sq^2==3 exactly: ", lift(Mod(sq,Pm)^2)==3);
ok=1; for(i=1,#V, my(E=expr(V[i]), Q=algdep(V[i],16)); if(lift(subst(Q,'x,Mod(E,Pm)))!=0, ok=0; print("inexact ",i))); print("all Ex exact: ",ok);
