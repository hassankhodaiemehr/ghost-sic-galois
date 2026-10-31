# Problem 1 (SICs) - session log, 2026-10-02

## A. Complexified SIC variety: complete complex solution counts (ghost_count.jl)
Unknowns x, y in P^{d-1} (y replaces conj(x)); equations (d+1)(y^T D_p x)(y^T D_p^dag x) = (y^T x)^2, p != 0.
All isolated solutions via total-degree/polyhedral homotopy (HomotopyContinuation.jl), filtered against the full system.

| d | paths | complex solutions | SIC fiducials (y = conj x) | non-real ("ghost") |
|---|-------|-------------------|----------------------------|--------------------|
| 4 | 560 | 512 | 256 | 256 |
| 5 | 14,614 | 4,000 | 2,000 | 2,000 |
| 6 | 218,343 | 6,912 | 3,456 | 3,456 |

Observation (conjecture G): V_d(C) = SICs ⊔ ghosts, with equal numbers; no other complex solutions.
If true for all d, then SIC existence <= [V_d(C) nonempty] + [reality lemma].
V_d(C) nonempty <= V_d has points over F̄_p for infinitely many p (Nullstellensatz + reduction mod p).

## B. Counting barrier for "soft" proofs
Restricting to a symmetry G scales unknowns (~2d/|G|) and equations (~d^2/|G|) equally, so the system stays overdetermined by ~d/2.
So intermediate-value, degree, or fixed-point arguments on symmetric slices can't work for large d. A proof must use non-generic algebraic structure.

## C. Heisenberg groups of non-cyclic abelian groups, d = 16..36 (hsearch.py, big.py; big.log)
L-BFGS on the frame potential, 60 random starts per group. Controls: Z_d is found for d <= 20, 25, 27; Hoggar (Z2^3) is found 30/30.
No SIC was found for any non-cyclic group: (4,4),(2,8),(2,2,4),(2^4),(3,6),(2,10),(2,12),(2,2,6),(5,5),(3,9),(3,3,3),(2,14),(2,16),(4,8),(2,2,8),(2,4,4),(2,2,2,4),(2^5),(2,18),(3,12),(6,6),(2,2,9),(2,3,6).
The best gaps are 0.026-0.27. Evidence is weak for d >= 24, where the cyclic control also failed at 60 starts.
