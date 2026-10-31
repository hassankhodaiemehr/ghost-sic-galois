Certified reciprocity law and the inert prime-power family (paper Secs. VI.G-H, Lemma VI.43, Thm. VI.46).

1. Ghost search (numerical, double precision):
     python3 findghost_z.py d seed      -> ghost{d}.txt   (x restricted to an eigenspace of the Clifford
                                                         unitary of F=[[0,-1],[1,-1]]; U_P-ghost, non-Hermitian)
2. Newton refinement to DIG digits (PARI/GP):
     ./run_ref.sh d DIG                 -> eps{d}.txt     ([labels, normalized overlaps eps_p])
   (uses grefine.gp and refine_drv.gp)
3. Exact field and reciprocity test:
     echo 'd=7; DIG=1500; QMAX=3000; NDEG=6; read("recip.gp")' | gp -q -s 2G > rec7.log
   NDEG = |Cl(d inf_-)| (expected [M_eps:F]).  recip.gp reads cert.gp at the end.
   Steps: stabilizer of (e_p) in GL2(Z/d); algdep + polredbest; every e_p as an exact element of
   K = Q(e_(0,1)) via nfroots; exact checks (e_p e_-p = 1, sqrt(D0) in K, Gal(K/F) abelian, kappa in
   Gal(K/F)); Frobenius at every prime ideal above q <= QMAX versus e_p -> e_{phi(varpi) p} for every
   generator varpi positive at the chosen place, for both order-3 stabilizer elements (Lam, Lam^-1) and
   both places (inf+, inf-); rnfconductor of K/F; ray classes hit by the tested primes.
   Exact expressions are saved in expr{d}.gp as [defining polynomial, Ex, theta].
   Data used in the paper: eps5.txt (800 digits), eps7.txt (1500), eps11.txt (6000), eps15.txt (4000).
   Logs: rec5.log, rec7.log, rec11.log, rec15.log (rec15_full.log contains the per-prime lines of the
   failing alternative conventions).
4. Family P:
     gp -q fam.gp  -> fam.log, family.txt  (members d<=2000: inertness, (C_d), class numbers, 2-rank)
     gp -q ker.gp  -> ker.log              (involutions in the image of P^m/P_{m,1} in Cl(d inf_-), d<=149)
