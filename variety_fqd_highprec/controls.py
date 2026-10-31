from hsearch import search
import time
for shape in [(7,),(2,2),(2,2,2),(2,4),(8,),(3,3),(9,),(2,6),(12,)]:
    t=time.time(); H,gap,hits,_=search(shape,starts=30)
    print(shape,'d=',H.d,'best gap=%.2e'%gap,'hits',hits,'/30 %.1fs'%(time.time()-t),flush=True)
