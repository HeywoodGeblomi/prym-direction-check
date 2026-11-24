# PREREG: fresh multi-twist closed-geodesic survey on the D=5 Prym eigenform in H(6)

Written before any word enumeration or spectrum computation. Construction checks (surface.py, checks.py,
rm_check.py, symmetries.py) were run before this file; they touch no word data.

## Surface (fixed, from the construction step)
- Thurston-Veech surface on a 10-rectangle grid: r=[7,9,1,8,6,4,5,0,3,2], u=[1,5,6,4,3,0,8,9,2,7],
  Prym involution sigma=[8,6,5,7,9,2,1,3,0,4]. H(6) (one zero of order 6), involution has exactly 2 fixed points.
- 4 horizontal and 4 vertical maximal cylinders, intersection matrix E (det -1), mu = phi^2, all moduli equal mu.
- Homology basis: the 8 cylinder cores (Z-basis since det E = -1). H1+ = sigma-invariant part (rank 4, saturated
  basis c_h0-c_h2, c_h1-c_h3, c_v1-c_v0, c_v3-c_v2 as computed). Real multiplication on H1- verified to be exactly O_5.
- Generators: H = horizontal Thurston-Veech multitwist (derivative [[1,mu],[0,1]]), V = vertical multitwist
  (derivative [[1,0],[-mu,1]]). Normalized by diag(mu,1): H=[[1,phi^4],[0,1]], V=[[1,0],[-1,1]], entries in Z[phi]=O_5.
  The group <H,V> is free (ping-pong, |mu^2| > 4). It may be a proper subgroup of the full Veech group; this is stated
  as a limit, not tested.

## Enumeration (primary run, the only run)
- Alphabet: letters H^k and V^k with k in {-4,-3,-2,-1,1,2,3,4} (negative shears included, so inverses and
  reversals are reachable).
- Words: all reduced alternating words (no two adjacent letters of the same type) of length 1 to 4, starting with
  either type. Count = 2*(8+64+512+4096) = 9360 directed words. Order of multiplication: word [l1,...,ln] means the
  affine map l1 o l2 o ... o ln (matrix product left to right).
- Hyperbolic: |trace(derivative)| > 2, computed exactly in Z[phi].

## Exact 2:1 criterion (raw H1+ spectrum)
- M_plus = action of the word on the saturated integral basis of H1+ (4x4 integer, symplectic for the restricted form).
- Let p(t) = charpoly(M_plus). Positive exponents l1 >= l2 >= 0 are log of the two largest root moduli.
- Exact 2:1 means l1 = 2*l2 > 0. Screen numerically at 60 digits (|l1/l2 - 2| < 1e-40), then confirm exactly:
  some root of p of maximal modulus is also a root of q(t) or q(-t), where q = charpoly(M_plus^2) (so mu1 = +-mu2^2).
  Hits that pass the numeric screen but not the exact identity are reported separately, not counted.
- Raw spectrum reported as [1, l2/l1, -l2/l1, -1]; exponents along the geodesic reported as l/log(rho(A)).

## Dedupe rules (applied to exact hits)
1. Directed hits: number of enumerated words that are exact 2:1.
2. Up to cyclic rotation: cyclically reduce each word (merge first and last letters while same type, drop zero
   exponents), then take the lexicographically least rotation. Since <H,V> is free this is exactly conjugacy in <H,V>.
3. Up to rotation and reversal: also identify w with w^-1 (reverse order, negate exponents), canonical form = least
   over all rotations of w and of w^-1.
4. Up to surface symmetry, reported separately: the only symmetries found are translation identity, the Prym
   involution (derivative -I, central, acts trivially on classes), and two orientation-reversing reflections
   (derivatives [[0,1],[1,0]] and its negative) which conjugate H^k -> V^-k, V^k -> H^-k. Class under rule 3 plus this
   letter map. These reflections send a geodesic to its mirror image, so this count is a different notion than
   "same closed geodesic"; it is reported only as a separate column.
- Consistency checks on every merged pair: equal |trace(A)| and equal charpoly of the full 8x8 action
  (conjugate maps must agree; a mismatch is a bug and is reported).
- Also report counts of distinct traces and distinct 8x8 charpolys among hits (necessary-condition proxies).

## Verdict rule
- Report whether exact 2:1 occurs at all on this D=5 surface in this word range, and the counts under 1 to 4.
- "20 stays 20" or "becomes 10" would only be meaningful if this run reproduced 20 directed hits; otherwise report the
  numbers as they come out. The original survey's 20 cannot be recounted (its word list is gone).

## Recorded example test (separate from the survey)
- Kenyon-Smillie / McMullen: trace field = Q(tr A) for any hyperbolic A in a Veech group, and it is Q only for torus
  covers. Report what trace 52 implies. Also try simple readings of the word [(1,3),(2,4),(2,4),(3,2),(0,2)] in <H,V>.
