# RESULTS: where the recorded example (A, M_plus) actually lives

Prereg: PREREG.md, sha256 ff6cab2a6edeabc8018b35cfad80f32b15491a1614eb9821e4544532137ccb69, written 2026-11-23 14:55:42 EST,
before the survey ran. Survey (survey.py) run once, all seven configurations.

## 1. Home of the example
- Exact reproduction (toy_model.py): take a 4-dim lattice c0,c1 (horizontal-type cores), c2,c3 (vertical-type cores)
  with form Om(c0,c2)=Om(c1,c3)=-1 and nothing else. Letter (i,k) = k-fold transvection on core i. Derivative = A_h(k)
  for i in {0,1}, A_v(k) for i in {2,3}, first letter acts first. The word [(1,3),(2,4),(2,4),(3,2),(0,2)] then gives
  exactly A=[[31,10],[65,21]] (= A_h(2) A_v(10) A_h(3)) and exactly the recorded M_plus.
- What that lattice is: two separate hyperbolic planes, i.e. H1 of two disjoint square tori. M_plus splits into
  block P=[[-15,-2],[8,1]] on (c0,c2), trace -14, and block Q=[[1,-3],[2,-5]] on (c1,c3), trace -4.
- Two problems make it formal, not a translation-surface map:
  a. Each letter twists one cylinder only. A real affine multitwist with derivative A_h(3) twists every cylinder in
     that direction, but letter (1,3) leaves block P untouched.
  b. Sign mismatch: the transvections are all the same handedness, while A_h(k) and A_v(k) with k>0 are twists of
     opposite handedness. On a single torus the homology action equals the derivative; the derivative-consistent blocks
     would be A_h(2)A_v(8) (trace 18) and A_v(2)A_h(3) (trace 8), whose exponent ratio is 1.399, not 2.
- Genuine search (home_search.py, home_search_out.json): 105 Prym-type origamis in H(6), 8 to 16 squares, words in the
  true horizontal and vertical multitwists (lengths 2, 4, 6). Words with |trace| 52 occurred on 11 of them; none had
  H1+ char poly (t^2+14t+1)(t^2+4t+1). Seen instead: (t-1)^2(t^2-4t+1), (t-1)^2(t^2-20t+1), (t-1)^2(t^2+28t+1),
  (t^2+16t+1)(t^2+34t+1), (t^2+16t+1)^2. No genuine surface found; the smallest honest "home" is the two-torus model.
- Reconciling the Q(sqrt3): any 2x2 integer block with eigenvalues in Q(sqrt3) has eigenvalue a power of the unit
  2+sqrt3 (unit group of rank 1). Traces 4, 14, 52 are the first three: (2+sqrt3)^1, ^2, ^3. Here 4 = |2-6|, 14 = |2-16|,
  52 = 2+50, small-number coincidences of the k values. In a real Prym surface H1+ has no tautological part, so A's
  eigenvalue does not have to appear there at all; the shared unit is coincidence, not a link.

## 2. Survey on the model (survey_out.json)
| Config | words | hyperbolic A | M_plus hyperbolic | two positive | band 1.9-2.1 | exact 2:1 directed | up to rotation | +reversal | up to symmetry | distinct M_plus |
|---|---|---|---|---|---|---|---|---|---|---|
| C1 len 5, k 1..4 | 1048576 | 983040 | 821860 | 109440 | 6720 | 4320 | 864 | 864 | 216 | 216 |
| C2 len 5, k 2..4 | 248832 | 233280 | 216420 | 51840 | 2400 | 960 | 192 | 192 | 48 | 80 |
| C3 len 1..5, k 1..4 | 1118480 | 1043584 | 862760 | 110976 | 6816 | 4416 | 888 | 888 | 222 | 216 |
| C4 len 1..4, k 1..4 | 69904 | 60544 | 40900 | 1536 | 96 | 96 | 24 | 24 | 6 | 16 |
| C5 random 10k, len 5, k 1..4 | 10000 | 9360 | 7766 | 967 | 68 | 43 | 40 | 40 | 35 | 35 |
| C6 random 10k, len 5, k 2..4 | 10000 | 9389 | 8714 | 2097 | 96 | 39 | 33 | 33 | 26 | 30 |
| C7 len 4, k +-1..+-4 | 1048576 | 701824 | 608432 | 55296 | 2496 | 2496 | 624 | 312 | 78 | 416 |
| Original survey | ? | 8095 | 5504 | 850 | 35 | 20 (46 matrices) | | | | |
- Merge consistency violations: 0 in every configuration.
- Reproduced: no. No configuration gives 8095. Closest by rates is C5 (random 10k): two-positive share 10.3% vs 10.5%,
  but M_plus-hyperbolic share 83% vs 68%, band 0.73% vs 0.43%, exact 0.46% vs 0.25%.
- Only C7 can reach inverses; there reversal halves the count (624 -> 312), the 20 -> 10 effect, because every
  hit's inverse is also a hit.

## 3. Is exact 2:1 here interesting?
- Exponents along one closed geodesic are log eigenvalue ratios of one matrix. They are not the Kontsevich-Zorich
  exponents of the Teichmuller curve, which are averages (EKZ sum formula for square-tiled surfaces, arXiv:1112.5872).
- On a square-tiled surface the affine group acts on H1 by integer matrices; H1 splits into the tautological plane
  and the rest, often further by symmetry (Matheus-Yoccoz-Zmiaikou, arXiv:1207.2423). Whenever two 2x2 integer blocks
  have eigenvalues in the same real quadratic field, the exponent ratio is forced to be rational, because that field's
  units are powers of one fundamental unit. 2:1 is just the smallest case (traces t and t^2-2).
- In the model it is common: 4320 directed hits per 1M words in C1.
- Verdict: expected. Not a finding about D=5, D=12 or the 4/7, 2/7 split (which is known: Yu-Zuo plus Moller, as the
  repo itself says).

## Harness line
Held up: the recorded A and M_plus are reproduced exactly by a formal two-torus transvection model with a sign
mismatch, no genuine Prym origami among 105 tried reproduces them, the old 8095/850/35/20 survey was not reproduced,
and exact 2:1 in this setting is forced by unit arithmetic, so it is expected, not a big finding.
