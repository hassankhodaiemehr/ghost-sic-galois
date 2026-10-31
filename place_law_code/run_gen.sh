#!/bin/bash
k=$1
python3 refx.py $k 700 >/dev/null 2>&1
python3 - $k <<'PY' > inx_$k.gp
import sys
k=sys.argv[1]
rows=[l.split() for l in open(f'x_{k}.txt')]
print("L=["+",".join(f"[{r[0]},{r[1]},{r[2]},{r[3]},{r[4]}]" for r in rows)+"];")
PY
gp -q -s 2G <<GP
default(realprecision,680);
read("hs.gp"); read("gen.gp"); read("inx_$k.gp");
genanalyze(L);
GP
