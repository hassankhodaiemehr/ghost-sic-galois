# ghost-sic-galois

Code, data and logs for the paper

> **H. Khodaiemehr**, *From Ghost SICs to SIC-POVMs: A Galois Reduction of Zauner's Conjecture* (2026), arXiv: (to be added).

The latest version of the paper is in [`paper/sic_ghosts.pdf`](paper/sic_ghosts.pdf). The code and data are also distributed as ancillary files (`anc/`) of the arXiv submission. Each code directory has its own `README.txt`.

| Directory | Contents (paper section) |
|---|---|
| `variety_fqd_highprec/` | Complete solution counts of the complexified SIC variety for d = 4, 5, 6 (Julia, HomotopyContinuation.jl); parity test; abelian-Heisenberg search; finite quantum dilogarithm solvers; PSLQ / high-precision checks; Problem 5 check. Summary in `RESULTS.md`. |
| `place_law_code/` | Place-law and single-unit checks; class-field test for d = 4, 5 (`d5/`). |
| `rayclass_code/` | Ray-class tests for d = 9, 13 (d = 17 attempt); unique-involution census. |
| `reciprocity_code_d5_d9/` | Explicit reciprocity law, first exact tests (d = 5, 9; Sec. VI.G, Table VII). |
| `ramification_coverage/` | Condition (C_d) for 4 ≤ d ≤ 20000; coverage figure (Fig. 1); reduced FQD system. |
| `recip_v4/` | Certified reciprocity law for d = 5, 7, 11, 15 (Lemma VI.43) and the prime-power family P (Sec. VI.H, Theorem VI.46). |

## Requirements

- PARI/GP 2.15 or later
- Python 3 with numpy, scipy and mpmath
- Julia 1.9 or later with HomotopyContinuation.jl

File paths inside the scripts are relative to the directory that contains the script. Typical runs are documented in each `README.txt`. For example:

```bash
cd recip_v4
./run_ref.sh 7 1500                     # refine the d=7 ghost to 1500 digits -> eps7.txt
echo 'd=7; DIG=1500; QMAX=3000; NDEG=6; read("recip.gp")' | gp -q -s 2G > rec7.log
```

## Status of the mathematics

The paper does **not** prove Zauner's conjecture. It reduces SIC existence to an archimedean Stark-type condition. On an explicit infinite family it further reduces that condition to an unramified condition. Section VII.C of the paper lists what remains.

## License

The code is released under the MIT License (see `LICENSE`). Data files and logs may be reused with citation of the paper.
