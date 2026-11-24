# PREREG: word survey on the model that reproduces the recorded (A, M_plus)

Written before any survey enumeration. Done before this file: toy_model.py (reproduces the single recorded word),
home_search.py (looked for a genuine Prym H(6) origami home; none found). No survey counts have been computed.

## Model ("toy", the only candidate that reproduces both A and M_plus exactly)
- 4-dim lattice with basis c0,c1 (horizontal-type cores), c2,c3 (vertical-type cores), symplectic form
  Om(c0,c2)=Om(c1,c3)=-1, all other pairs 0 (two orthogonal hyperbolic planes, i.e. H1 of two separate tori).
- Letter (i,k): homology = transvection x -> x + k*Om(c_i,x)*c_i. Derivative = A_h(k)=[[1,0],[k,1]] if i in {0,1},
  A_v(k)=[[1,k],[0,1]] if i in {2,3}. A word acts letter by letter, first letter first.
- M_plus splits into block P on (c0,c2) and block Q on (c1,c3), each in SL(2,Z).

## Definitions (per word)
- hyperbolic_A: |tr A| > 2.
- hyperbolic_M_plus: at least one block has |trace| > 2.
- two_positive_exponents: both blocks have |trace| > 2.
- exponents l(X) = log of the larger eigenvalue modulus of block X = arccosh(|tr X|/2).
- exact 2:1 (exact integer test): both blocks hyperbolic and |tr big| = (tr small)^2 - 2, where big is the block with the
  larger |trace|. (This is exactly lambda_big = +-lambda_small^2.)
- band: 1.9 <= l_big/l_small <= 2.1 (both hyperbolic).
- n_matrices_21: number of distinct M_plus matrices among exact hits.

## Configurations (all run once, in this order)
- C1: all words of length exactly 5, i in {0,1,2,3}, k in {1,2,3,4}: 16^5 = 1,048,576 words.
- C2: all words of length exactly 5, k in {2,3,4}: 12^5 = 248,832 words.
- C3: all words of length 1 to 5, k in {1,2,3,4}: 1,118,480 words.
- C4: all words of length 1 to 4, k in {1,2,3,4}: 69,904 words.
- C5: 10,000 random words of length 5, k in {1,2,3,4}, Python random.seed(0), letters drawn uniformly i.i.d.
- C6: 10,000 random words of length 5, k in {2,3,4}, random.seed(0).
- C7: C1 alphabet but k in {-4..-1,1..4}, length exactly 4 (32^4 = 1,048,576), so inverses are reachable.

## Dedupe of exact hits (each configuration)
1. directed: word count.
2. up to rotation: least cyclic rotation of the word.
3. up to rotation and reversal: least over rotations of w and of w^-1 (reverse order, negate k). With positive-only
   alphabets w^-1 is never in the set, so 2 and 3 can only differ in C7.
4. up to symmetry: also the letter relabelings that preserve the model: block swap (0<->1, 2<->3) and type swap
   (0<->2, 1<->3), and their composition.
- Every merged class is checked for equal M_plus char poly and equal |tr A| (necessary for conjugacy).

## Reproduction rule
- "Reproduced" only if some configuration gives exactly hyperbolic_A = 8095 together with 850 two-positive and
  20 exact hits (or 46 matrices). Otherwise report the closest configuration by rates
  (5504/8095 = 0.680, 850/8095 = 0.105, 35/8095 = 0.0043, 20/8095 = 0.0025) and say not reproduced.
