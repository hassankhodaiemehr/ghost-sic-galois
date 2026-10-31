#!/bin/bash
d=$1; f=$2; prec=$3
python3 - $f <<'PY' > in5.gp
import sys
rows=[l.split() for l in open(sys.argv[1]) if not l.startswith('U_P')]
print("E=["+",".join(f"[{r[0]},{r[1]},{r[2]}]" for r in rows)+"];")
PY
gp -q -s 2G <<GP
default(realprecision,$prec);
read("an5.gp"); read("in5.gp");
L=analyzeGhost(E,$d);
\\gen
GP
