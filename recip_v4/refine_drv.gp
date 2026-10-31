\\ usage: gp -q -s 4G with d, DIG set: refine ghost{d}.txt and write eps{d}.gpv
read("grefine.gp");
default(realprecision, DIG+50);
{
my(f=readstr(Str("ghost",d,".txt")), re=eval(Str("[",strjoin(strsplit(f[1]," "),","),"]")), im=eval(Str("[",strjoin(strsplit(f[2]," "),","),"]")));
my(ov=ghostrefine(d,re,im,DIG), ps=List());
for(a=0,d-1, for(b=0,d-1, if(a||b, listput(ps,[a,b]))));
my(maximag=vecmax(apply(z->abs(imag(z)),ov))); print("max |imag eps| = ", maximag*1.);
write(Str("eps",d,".txt"), [Vec(ps), apply(real,ov)]);
print("saved ",#ov);
}
