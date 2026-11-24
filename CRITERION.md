# Direction-check criterion (written before reading any orbit data or dedupe results)

Date: 2026-11-23, written before fetching any data/*.json or docs from HeywoodGeblomi/prym-omega-ed6-h1plus.

Question: does the recorded count of "20 exact 2:1 closed geodesics" double-count geodesics?

Definition (gold standard). Two recorded closed geodesics are the SAME unoriented closed geodesic iff
their shear words are equal up to (a) cyclic rotation and/or (b) inversion (reverse the word and invert
each letter; A_h(k)^-1 = A_h(-k), A_v(k)^-1 = A_v(-k)). Equivalently their derivative matrices are
conjugate in the Veech / affine group, up to inversion (A ~ A^-1).

Counts to report, if per-orbit data exist:
  N_dir   = number of directed recorded orbits (raw rows)
  N_rot   = distinct up to cyclic rotation
  N_unor  = distinct up to cyclic rotation + reversal/inversion
Optional (c): also identify under the surface involution / affine symmetries only if the repo gives them.

Proxies, in decreasing strength, and I will say which I could apply:
  P1 words available: canonicalize each word by min over all rotations of word and of inverse word.
  P2 only matrices: A ~ B if B equals A, A^-1, or a conjugate by a known group element; weaker proxy:
     equal trace and (B == A or B == A^-1) exactly.
  P3 only traces / exponents: trace equality is NECESSARY but not sufficient for sameness, so it can give
     only an upper bound on collapse, never a proof of duplication.
If no per-orbit list exists, the verdict is "undetermined" and I name exactly what data are missing.

Verdict options: 20 stays 20; becomes 10 (exact pairing by reversal); some other number; or undetermined.
