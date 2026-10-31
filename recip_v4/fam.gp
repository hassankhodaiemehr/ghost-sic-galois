\\ family P = {p^k : p = 5 mod 6 prime}, d <= B : class number of F, ramification at inf-, 2-Sylow of Cl(d inf-) and of P/P_{m,1}-image
B=2000; L=List();
{forprime(p=5,B, if(p%6!=5,next); for(k=1,30, my(d=p^k); if(d>B,break); listput(L,d)))};
L=vecsort(Vec(L));
{res=List(); for(i=1,#L, my(d=L[i], Dl=(d-3)*(d+1), D0=core(Dl), bnf=bnfinit(y^2-D0,1), h=bnf.no,
   rts=bnf.roots, im=if(rts[1]<0,1,2), arch=vector(2,j,j==im), c1=bnrclassno(bnf,[d,arch]), c0=bnrclassno(bnf,d),
   inert=#idealprimedec(bnf,factor(d)[1,1])==1, cyc=bnrinit(bnf,[d,arch]).cyc, t2=#[c|c<-cyc, c%2==0]);
   listput(res,[d,D0,h,inert,c1==2*c0,t2]))};
print("members d<=",B,": ",#res);
print("all inert: ",vecmin(apply(r->r[4],Vec(res))), "   inf- ramified in all: ", vecmin(apply(r->r[5],Vec(res))));
print("h odd: ",#[r|r<-Vec(res),r[3]%2==1], "   list: ",[r[1]|r<-Vec(res),r[3]%2==1]);
print("2-rank of Cl(d inf-) = 1 (cyclic 2-Sylow): ",#[r|r<-Vec(res),r[6]==1]);
print("first 30: ", apply(r->[r[1],r[2],r[3]],Vec(res)[1..30]));
write("family.txt", res);
