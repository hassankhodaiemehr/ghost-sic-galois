#!/bin/bash
# usage: an.sh k maxd
k=$1; md=${2:-40}
python3 - $k <<'PY' > in_$k.gp
import sys
k=sys.argv[1]
rows=[l.split() for l in open(f't_{k}.txt')]
print("L=["+",".join(f"[{r[0]},{r[1]},{r[2]},{r[3]},{r[4]}]" for r in rows)+"];")
PY
timeout 1800 gp -q -s 2G <<GP
default(realprecision,850);
read("hs2.gp"); read("in_$k.gp");
summ(L,$md);
GP
