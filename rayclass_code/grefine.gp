ghostrefine(d, re, im, digits)={
  my(w=exp(2*Pi*I/d), tau=-exp(Pi*I/d), X=matrix(d,d,i,j,if((i-j-1)%d==0,1,0)), Z=matdiagonal(vector(d,k,w^(k-1))),
     UP=matrix(d,d,i,j,if((i-1+j-1)%d==0,1,0)), Ds=List(), DH=List(), x, y, z, it);
  for(a=0,d-1, for(b=0,d-1, if(a||b, my(D=tau^(a*b)*X^a*Z^b); listput(Ds,D); listput(DH,conj(D)~))));
  Ds=Vec(Ds); DH=Vec(DH);
  x=vector(d,k,re[k]+I*im[k])~; x=x/x[1]; y=conj(UP*x); my(s=y~*x); y=y/s;
  z=concat(x,y);
  for(it=1,60,
    my(x=z[1..d], y=z[d+1..2*d], s=y~*x, F=vector(#Ds+2), J=matrix(#Ds+2,2*d));
    for(p=1,#Ds, my(A=Ds[p],B=DH[p], Ax=A*x, Bx=B*x, yA=y~*A, yB=y~*B, a=y~*Ax, b=y~*Bx);
      F[p]=(d+1)*a*b-s^2;
      for(k=1,d, J[p,k]=(d+1)*(yA[k]*b+a*yB[k])-2*s*y[k]; J[p,d+k]=(d+1)*(Ax[k]*b+a*Bx[k])-2*s*x[k]));
    F[#Ds+1]=x[1]-1; J[#Ds+1,1]=1;
    F[#Ds+2]=s-1; for(k=1,d, J[#Ds+2,k]=y[k]; J[#Ds+2,d+k]=x[k]);
    if(normlp(F,oo)<10^(-digits+15), break);
    my(JH=conj(J)~); z=z+matsolve(JH*J, -(JH*F~)));
  my(x=z[1..d], y=z[d+1..2*d], s=y~*x);
  vector(#Ds,p, sqrt(d+1)*(y~*Ds[p]*x)/s)
}
