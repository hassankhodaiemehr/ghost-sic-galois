Ray-class tests of the overlap-unit field, d = 9, 13 (d = 17 pending). PARI/GP 2.15+, Python 3 (numpy, scipy).

findghost.py, findsic.py   numerical U_P-ghost / WH SIC by least squares (inputs ghost*.txt, sic*.txt)
grefine.gp, refine.gp      Newton refinement on V_d (ghost / SIC) to high precision; gin*.gp, in*.gp drivers
g9.gp, g13x.gp, g17x.gp    refined ghosts: d=9 (2000 digits), d=13 and d=17 (20500 digits)
e9.gp, e9h.gp              refined WH SIC for d=9 (SIC-side mirror check)
d9*.gp, ganal.gp, an.gp    d=9: minimal polynomials, rnfconductor, Galois group over Q(sqrt15)
d13a-e.gp, d13polys.gp     d=13: trace polynomial Q13 (degree 96), conductor (13) of the trace field
d13Rt.gp, d13f.gp          d=13: Frobenius test (Frobenius order vs ray class order in Cl((13)inf-), kappa test)
d17.gp, d17f.gp            d=17 attempt (not identified at 20500 digits)
sl9.gp, sl13.gp            place-law sign check on all conjugates (uses hs2.gp: signlaw)
rc.gp                      ray class groups Cl(d'inf-), Cl(d') for d = 4,5,9,13,17
twosyl.gp                  2-Sylow of Cl(d'inf-) for 4<=d<=150 (unique-involution list)
twosyl2.gp                 census of cyclic 2-Sylow with inf- ramified, 4<=d<=2000 (24 dimensions)
comp.gp                    helper (compositum experiments; unused in the paper)
*.log                      run logs
