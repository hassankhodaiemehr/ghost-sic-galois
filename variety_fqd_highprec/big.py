from hsearch import search
import time, sys, json
shapes=[(16,),(4,4),(2,8),(2,2,4),(2,2,2,2),(18,),(3,6),(20,),(2,10),(24,),(2,12),(2,2,6),(25,),(5,5),(27,),(3,9),(3,3,3),(28,),(2,14),
        (32,),(2,16),(4,8),(2,2,8),(2,4,4),(2,2,2,4),(2,2,2,2,2),(36,),(2,18),(3,12),(6,6),(2,2,9),(2,3,6)]
out=[]
for shape in shapes:
    t=time.time(); H,gap,hits,_=search(shape,starts=int(sys.argv[1]))
    r=dict(shape=shape,d=H.d,gap=gap,hits=hits); out.append(r)
    print(shape,'d=',H.d,'best gap=%.2e'%gap,'hits',hits,'%.0fs'%(time.time()-t),flush=True)
json.dump(out,open('abelian_results.json','w'))
