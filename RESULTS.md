# Direction check results (after CRITERION.md, sha256 f47d93c51b3af8d2c9cf2ac70c973227b6976e2b347ac81cc30e992fe764f2a1)
Repo HeywoodGeblomi/prym-omega-ed6-h1plus @ c3f09f506f51 (read via API only; repo is public per API, not private).
Per-orbit list of the 20: NOT in repo, not in any of its 11 commits. No code in repo. data/x11_affine_maps.json is cited by paper/computational_note_d6.md but never committed.
Only 3 of the 20 are identifiable (sample_exact_21 rho_A 52.0, 44.0, 148.0 -> traces 52, 44, 148). Distinct traces => 3 pairwise distinct unoriented closed geodesics (trace is invariant under rotation, inversion, conjugation).
Counts: directed 20 (as recorded); up to rotation: undetermined (>=3); up to rotation+reversal: undetermined (>=3). Proxy applied: P3 trace only.
Example: A=[[31,10],[65,21]] tr 52, rho=(2+sqrt3)^3; M+ eigenvalues -(2+sqrt3)^2, -(2+sqrt3) and inverses -> ratio exactly 2:1, spectrum symmetric so A^-1 is also an exact 2:1 hit.
Naive shear product of example word (even index=h or =v) gives trace 148, not 52 -> word-to-A map not reconstructable from repo.
