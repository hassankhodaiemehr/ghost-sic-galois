default(realprecision,3050);
read("an.gp"); read("e9h.gp");
{
my(e=e9[10], z=exp(2*Pi*I*12/36)*e^2, t=real(z+1/z), P=algdep(t,36), f=factor(P)[,1], Q);
for(i=1,#f, if(abs(subst(f[i],'x,t))<10^-2000, Q=f[i]));
my(Pz=simplify('x^36*subst(Q,'x,'x+1/'x)));
print("irreducible ",polisirreducible(Pz));
print("t field: ",relcond(Q,t,15));
print("z field: ",relcond(Pz,z,15));
my(bnf=bnfinit(y^2-15,1)); print("Cl(9 * place2) = ",bnrinit(bnf,[9,[0,1]]).cyc, "  Cl(9) = ",bnrinit(bnf,9).cyc);
write("d9Pz.gp","Pz9=",Pz,";");
}
