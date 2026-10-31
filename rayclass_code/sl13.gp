\\ Place-law sign check for d=13 on the degree-96 trace polynomial Q13 (sqrt(35)).
read("hs2.gp"); read("d13polys.gp");
default(realprecision,600);
{ my(r=polroots(Q13), i0=0); for(i=1,#r, if(abs(imag(r[i]))<1e-100 && real(r[i])>2, i0=i; break));
print(signlaw(Q13,35,real(r[i0]))); }
