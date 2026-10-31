Explicit reciprocity law (Sec. VI.G, Conjecture VI.40, Table VII): exact tests for d=5 and d=9.
PARI/GP 2.15+. Run as:  gp -q -s 8G < script.gp
Paths inside the scripts are relative: run them from this folder (g9.gp is read from ../rayclass_code).
g9.gp (2000-digit d=9 ghost) and grefine.gp are in rayclass_code.

d=5
  g5.gpv      ghost overlaps [p1,p2,value], 400 digits
  all5.gp     stabilizer, exact expressions, Frobenius test at split primes q<=3000
              (key "211-" = phi(sqrt3)=Lambda^-1-2I, no transpose/inverse, positive at inf_-: 594/594)
  ver5.gp     exactness check of the expressions
  inert5.gp   inert primes: Frob = qI, 219/219

d=9
  gal9.gp     minimal polynomial of e_(1,0) (degree 72) and its 36 automorphisms -> aut9.gp
  all9h.gp    exact expressions of all 80 squared overlaps -> ex9_all.gp
  all9i.gp    Frobenius test, q<=4000 (log all9i.log): key "211-" 122/122 principal split,
              inert qI 276/276, 136 non-principal primes move the ghost to its partner
