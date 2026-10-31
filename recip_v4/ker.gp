L=[5,11,17,23,25,29,41,47,53,59,71,83,89,101,107,113,121,125,131,137,149];
inv2(cyc,v)=for(i=1,#cyc, if((2*v[i])%cyc[i],return(0))); 1;
{foreach(L,d, my(Dl=(d-3)*(d+1), D0=core(Dl), bnf=bnfinit(y^2-D0,1), rts=bnf.roots, im=if(rts[1]<0,1,2), arch=vector(2,j,j==im),
   bnr=bnrinit(bnf,[d,arch],1), cyc=bnr.cyc, gens=bnr.bid.gen, G=vector(#gens,j,bnrisprincipal(bnr,gens[j],0)));
   my(S=Map(), Q=List(), z=vector(#cyc)~, cnt2=0); mapput(S,z,1); listput(Q,z); my(i=1);
   while(i<=#Q, my(v=Q[i]); i++; for(j=1,#G, my(w=vector(#cyc,t,(v[t]+G[j][t])%cyc[t])~); if(!mapisdefined(S,w), mapput(S,w,1); listput(Q,w))));
   for(t=1,#Q, if(inv2(cyc,Q[t]), cnt2++));
   print(d,": h=",bnf.no,"  |Cl_m|=",vecprod(cyc),"  |image of P| = ",#Q," = (p^2-1)/3? ",#Q==(factor(d)[1,1]^2-1)*d^2/factor(d)[1,1]^2/3,"   #elements of order<=2 in image: ",cnt2))}
