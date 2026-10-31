python3 - <<'PY' > in5x.gp
rows=[l.split() for l in open('g5.txt')]
print("L=["+",".join(f"[0,50,0,{i+1},{r[2]}]" for i,r in enumerate(rows))+"];")
PY
gp -q -s 4G <<GP
default(realprecision,380);
read("hs2.gp"); read("gen.gp"); read("in5x.gp");
L=[L[1],L[2],L[3],L[4],L[5],L[6],L[7],L[8],L[9],L[10],L[11],L[12],L[13],L[14],L[15],L[16],L[17],L[18],L[19],L[20],L[21],L[22],L[23],L[24]];
P0=findpol(L[1][5],64); print("deg x_u0 = ",poldegree(P0));
y=Mod('x,P0); ok=1;
for(j=2,#L, if(abs(abs(L[j][5])-1)<1e-50,next); Pj=findpol(L[j][5],64); f=expr(P0,L[1][5],Pj,L[j][5]); if(f==0, print("v=",j," not in Q(x_u0)"); ok=0; next); if(subst(f,'x,y)*subst(f,'x,y^-1)!=1, print("v=",j," kappa fails"); ok=0));
print("G1&G2 for u0=(0,1): ",ok);
GP
