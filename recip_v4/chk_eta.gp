default(realprecision,300);
idx(P,d,p)=my(a=p[1]%d,b=p[2]%d); for(i=1,#P, if(P[i][1]==a&&P[i][2]==b, return(i))); 0;
chk(d,eta,Lam)=my(W=read(Str("eps",d,".txt")), P=W[1], V=apply(z->z^2,W[2]), Dl=(d-3)*(d+1), D0=core(Dl), k=sqrtint(Dl/D0), c=Vec(lift(eta)), X=(2*Lam^-1-(d-1))/k, G, ok=1); if(#c==1,c=[0,c[1]]); G=Mod(c[2],d)*matid(2)+Mod(c[1],d)*X; for(i=1,#P, if(abs(V[idx(P,d,lift(G*P[i]~))]-V[i])>1e-200, ok=0)); print("d=",d," eta=",lift(eta)," phi(eta)=",lift(G)," det=",lift(matdet(G))," in stabilizer: ",ok);
chk(7, Mod(-1-y,y^2-2), Mod([3,5;3,3],7));
chk(15, Mod(2+y,y^2-3), Mod([0,14;1,14],15));
chk(7, Mod(-1-y,y^2-2), Mod([3,2;4,3],7));
