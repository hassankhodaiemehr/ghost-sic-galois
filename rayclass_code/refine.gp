\\ d, file sic<d>.txt -> refined psi; uses D_p = tau^(p1 p2) X^p1 Z^p2
sicrefine(d, re, im, digits)={
  my(w=exp(2*Pi*I/d), tau=-exp(Pi*I/d), X=matrix(d,d,i,j,if((i-j-1)%d==0,1,0)), Z=matdiagonal(vector(d,k,w^(k-1))), Ds=List(), psi, it, n, o, g, J, del);
  for(a=0,d-1, for(b=0,d-1, if(a||b, listput(Ds, tau^(a*b)*X^a*Z^b))));
  Ds=Vec(Ds);
  psi=vector(d,k,re[k]+I*im[k])~;
  for(it=1,40,
    my(prec=getlocalprec());
    n=real(conj(psi)~*psi);
    g=vector(#Ds); J=matrix(#Ds,2*d-2);
    for(p=1,#Ds, my(Dp=Ds[p], Dpsi=Dp*psi, psiD=conj(psi)~*Dp); o=conj(psi)~*Dpsi;
       g[p]=(d+1)*norm(o)-n^2;
       for(k=2,d, my(da=Dpsi[k]+psiD[k], db=-I*Dpsi[k]+I*psiD[k]);
          J[p,2*(k-1)-1]=(d+1)*2*real(conj(o)*da)-2*n*2*real(psi[k]);
          J[p,2*(k-1)]  =(d+1)*2*real(conj(o)*db)-2*n*2*imag(psi[k])));
    my(res=normlp(g,oo)); if(res<10^(-digits+10), break);
    del=matsolve(J~*J, -(J~*g~));
    for(k=2,d, psi[k]+=del[2*(k-1)-1]+I*del[2*(k-1)]));
  n=real(conj(psi)~*psi);
  [psi, vector(#Ds,p,sqrt(d+1)*(conj(psi)~*Ds[p]*psi)/n), it]
}
