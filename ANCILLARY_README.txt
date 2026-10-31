Ancillary files for "From Ghost SICs to SIC-POVMs: A Galois Reduction of Zauner's Conjecture"
(H. Khodaiemehr), arXiv version v4.
The same files are maintained at https://github.com/hassankhodaiemehr/ghost-sic-galois

Software: PARI/GP 2.15+, Python 3 (numpy, scipy, mpmath), Julia 1.9+ (HomotopyContinuation.jl).
File paths in the scripts are relative to the folder containing the script.

variety_fqd_highprec/   Complete solution counts of V_d (Julia homotopy), parity test, abelian-Heisenberg
                        search, finite quantum dilogarithm solvers, PSLQ/high-precision checks, Problem 5.
                        Results summary: RESULTS.md.
place_law_code/         Place law and single-unit checks; d=4,5 class-field test (folder d5).
rayclass_code/          Ray-class tests d=9,13 (d=17 attempt), unique-involution census. See README.txt there.
reciprocity_code_d5_d9/ Explicit reciprocity law, first computation (d=5, d=9). See README.txt there.
ramification_coverage/  Condition (C_d) for 4<=d<=20000 (place_ramification.py + CSV), coverage figure
                        (plot_route_c_coverage.py + CSV), reduced FQD system (fqd_sic.py + CSV).
recip_v4/               New in v4 (Secs. VI.G-H): Zauner-symmetric ghost search, Newton refinement,
                        certified reciprocity law in d=5,7,11,15(,19), and the family P = {p^k : p = 5 mod 6}.
                        See recip_v4/README.txt.
