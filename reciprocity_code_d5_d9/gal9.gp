read("../rayclass_code/g9.gp");
default(realprecision,2000); TOL=1e-1700;
e=vector(80,i,real(g9[i])^2);
mp(v,n)={my(Q=algdep(v,n),f=factor(Q)[,1]); for(k=1,#f, if(abs(subst(f[k],'x,v))<TOL, return(f[k]))); 0};
Pm=mp(e[9],72);
t0=getabstime(); A=nfgaloisconj(Pm); print("#aut ",#A," ms ",getabstime()-t0);
write("aut9.gp","Pm9=",Pm,";\nAut9=",A,";");
