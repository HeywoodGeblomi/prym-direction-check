# RESULTS: fresh D=5 rebuild of the multi-twist closed-geodesic survey

Preregistration: PREREG.md, sha256 1113290fd225bafd7188b1642572a3d95b60b7fde806810812278b8132747937,
written 2026-11-23 14:33:20 EST, before any word was enumerated. Survey run once (run_survey.py), right after.

## Headline
- Exact 2:1 does not occur on the genuine D=5 surface in the preregistered range.
  Directed hits 0, up to rotation 0, up to rotation and reversal 0, up to symmetry 0.
- Hyperbolic words: 9216 of 9360 (original survey reported 8095, so the scale is comparable).
- So the rebuild does not give 20, and does not give 10. It gives 0.
- The recorded example cannot come from a D=5 surface, or from a D=12 one. Its trace 52 is a whole number,
  and by Kenyon-Smillie / McMullen (trace field = Q(tr A) for any hyperbolic A; field Q only for torus covers)
  a pseudo-Anosov with trace 52 lives on a square-tiled (torus-cover) surface.

## 1. The surface (built and checked)
- Construction: Thurston-Veech on a 10-rectangle grid (an origami deformed by Perron-Frobenius lengths).
  r=[7,9,1,8,6,4,5,0,3,2], u=[1,5,6,4,3,0,8,9,2,7]. Found by a random search over origamis built with a
  built-in involution (search2.py); choice of this grid is mine.
- Stratum H(6): one zero, cone angle 14 pi. Genus 4.
- Prym involution sigma=[8,6,5,7,9,2,1,3,0,4] (rotation by pi), exactly 2 fixed points, quotient genus 2.
- Horizontal and vertical directions each split into 4 maximal cylinders (every cylinder boundary hits the zero).
  Intersection matrix E=[[1,0,0,1],[1,1,0,1],[0,1,1,0],[1,1,1,0]], det E = -1.
- Moduli: E E^T has top eigenvalue (7+3 sqrt5)/2 = phi^4, so mu = phi^2 and every cylinder has modulus mu.
  Heights (1, phi, 1, phi), circumferences (phi^2, phi^3, phi^2, phi^3). After rescaling, the horizontal
  cylinders are two unit squares and two squares of side (sqrt5-1)/2, which is exactly Lanneau-Nguyen
  genus-four paper (arXiv:1802.04879) Prop 2.4, the D=5 surface.
- Generators: H = horizontal multitwist (derivative [[1,mu],[0,1]]), V = vertical multitwist ([[1,0],[-mu,1]]).
  Normalized by diag(mu,1): H=[[1,phi^4],[0,1]], V=[[1,0],[-1,1]], entries in Z[phi] = O_5. Passed.
- Homology: the 8 cylinder cores are a Z-basis of H1 (Gram det = (det E)^2 = 1). Intersection form
  [[0,E],[-E^T,0]]. Twist action T_h(c_v_j) = c_v_j + sum_i E_ij c_h_i, and the mirror image for V.
- Checks (checks_out.json, rm_out.txt): T_h, T_v symplectic: yes. Prym involution is an involution, symplectic,
  commutes with both twists: yes. Period map intertwines each twist with its derivative, and the involution
  with -I: yes. H1+ (involution-fixed part) rank 4, saturated, periods of omega vanish on it: yes.
- Real multiplication: integer endomorphisms of H1- that scale omega form exactly Z[T] with T^2 = T + 1,
  T self-adjoint. That is O_5, discriminant 5 (not 20 or 45). With Prym, H(6) and McMullen's definition this
  puts the surface in Omega E_5(6), which Lanneau-Nguyen show is a single connected locus.
- Post-run check (diagnostics_out.txt): on every one of the 9216 hyperbolic words the top exponent on H1- equals
  log rho(A) (max relative error 5e-15), and all 8x8 actions are symplectic. So the homology action matches the
  affine maps.
- Symmetries (symmetries_out.txt): no translation automorphisms besides the identity; the Prym involution
  (derivative -I, central); no 90 degree rotations; two orientation-reversing reflections with derivatives
  [[0,1],[1,0]] and its negative (they send H^k to V^-k and V^k to H^-k).

## 2. Survey numbers (survey_out.json)
| Quantity | Count |
|---|---|
| Directed words enumerated (lengths 1 to 4, exponents +-1..+-4) | 9360 |
| Hyperbolic | 9216 |
| Two positive H1+ exponents | 7080 |
| Ratio in [1.9, 2.1] (post-run diagnostic, not preregistered) | 856 |
| Exact 2:1, directed | 0 |
| Exact 2:1, up to rotation | 0 |
| Exact 2:1, up to rotation and reversal | 0 |
| Exact 2:1, up to symmetry | 0 |
| Passed the numeric screen but failed the exact check | 0 |
- Closest ratio to 2 seen: 1.99810 (post-run diagnostic). Ratio quantiles 5%, 50%, 95%: 1.947, 2.854, 7.122.

## 3. The recorded example
- A=[[31,10],[65,21]]: det 1, trace 52, eigenvalues 26 +- 15 sqrt3 = (2+sqrt3)^(+-3).
- Recorded M_plus: charpoly (t^2+14t+1)(t^2+4t+1), eigenvalues -(2+sqrt3)^(+-2) and -(2+sqrt3)^(+-1).
  So the 2:1 there is exact, and it is arithmetic: both are powers of the same unit 2+sqrt3.
- D=5: the trace field of the D=5 surface is Q(sqrt5). A hyperbolic Veech element with trace 52 would make the
  trace field Q. Not consistent.
- D=12: the sqrt3 in the eigenvalues looks like D=12, but the same theorem applies: trace 52 forces trace field Q,
  while D=12 needs Q(sqrt3). Not consistent either.
- What it is consistent with: a torus cover (square-tiled surface), or integer shear algebra with no surface.
- Readings of the word in this D=5 group give traces 2083+933 sqrt5 or -173-75 sqrt5, never 52.
  Multiplying integer shears directly gives trace 148, also not 52. The word-to-A map is not recoverable.

## Limits
- <H,V> is free and may be a proper subgroup of the Veech group. A larger alphabet (other parabolic directions,
  or roots of H and V if they exist) could reach geodesics this run cannot. Rotation dedupe is exact conjugacy only
  inside <H,V>.
- Only one surface on the D=5 curve was used; that is fine for conjugacy-invariant spectra but the alphabet depends
  on it.
- One bounded range was run. Zero hits here does not prove exact 2:1 never occurs on D=5.

## Harness line
Held up: the D=5 surface is real (H(6), Prym, O_5 exactly, homology matches the affine maps on all 9216 words).
On it, exact 2:1 never appeared, and the recorded 2:1 example belongs to a square-tiled world, not D=5 or D=12.
