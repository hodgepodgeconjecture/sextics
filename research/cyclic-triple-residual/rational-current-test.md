# Testing rational lifts of the formal current

Research checkpoint, 2026-11-27, by OpenAI Codex for
`hodgepodgeconjecture`. This is a bounded rational-current test, not a
certificate or an unrestricted obstruction.

## Candidate chosen from the constructive calculation

The [formal construction](constructive-current.md) provides a leading
radial coefficient \(f_3=C_3/48\). At the seed \(J=1/(6i)\), so the
leading coefficient of \(g=f/J\) is \(iC_3/8\).

In normalized endpoint coordinates
\(u=k_0/k_1,v=k_2/k_1,w=k_3/k_1^3\), define

\[
\mathcal M=\{u^3,u^2v^2/w,uv^3/w,uv^2,v^3,v^3w/u\},
\qquad
\mathcal L=\{1,\xi/k_1,a,b,c,d/k_1^2,e/k_1,f\}.
\]

Here \(a,b,c,d,e,f\) are the existing line coefficients; the last \(f\)
in \(\mathcal L\) is the cubic line coefficient, not the radial unknown.
The 48 Laurent monomials \(\mathcal M\mathcal L\) allow endpoint poles
which the earlier polynomial screen omitted. The six monomial types in
\(\mathcal M\) were selected from the actual leading current, not by
raising a degree bound. Their leading coefficients are respectively

\[
-5i/288,\quad 3/(4\kappa^3),\quad 3i/(4\kappa^3),\quad
3i/\kappa^4,\quad -12/\kappa^4,\quad 9i/(8\kappa^3).
\]

Their linear combination is checked symbolically to equal \(iC_3/8\).
All candidate functions have weight zero under the original weighted
scaling, so \(b_i=Jk_i g_i\) has the required weight. Set \(b_1=0\).

Two classes were tested on \(k_1=1\):

1. Radial: \(g_0=g_2=g_3\), with 48 constant coefficients.
2. Independent components: each \(g_i\) belongs to this span independently,
   with 144 constant coefficients in total.

These are specific candidate classes. Neither is a consequence of a
complete pole-bound theorem; in particular they omit arbitrary nonlinear
dependence on the line coefficients and other denominator divisors.

## Full incidence-field evaluation

The tests solve all seven exact incidence equations at each p-adic
endpoint sample and differentiate those equations implicitly. They
evaluate the exact scalar source \(B\), not its leading homogeneous
term or a fitted single-ray approximation. The source deformation remains
\(F-sA\), evaluated at \(s=0\), with the fixed target unchanged.

For a field \(V\) equal to the local radial field (divergence 3), or
one component \(k_i\partial_i\) (divergence 1), the column for \(g\)
after clearing by \(18\Delta\xi/J\) is

\[
18\Delta\xi V(g)+
 [18\operatorname{div}(V)\Delta\xi-90\Delta V(\xi)-3\xi V(\Delta)]g.
\]

The right side is \(-2\Delta k_0^3/\xi^2-\xi\Delta_s\).
Both formulas are checked symbolically. Error precision is propagated
through rational evaluations and row reduction. All sample denominators
are nonzero on their Hensel branches.

| Class | Prime | Coefficient rank | Augmented rank | Final nonzero RHS valuation | Absolute precision |
| --- | ---: | ---: | ---: | ---: | ---: |
| Radial | 101 | 48 | 49 | 8 | 49 |
| Radial | 109 | 48 | 49 | 8 | 49 |
| Independent components | 101 | 144 | 145 | 9 | 47 |
| Independent components | 109 | 144 | 145 | 8 | 47 |

The respective normalized nonzero RHS residues are 16, 95, 66, and 57.
The coefficient ranks reach the total numbers of columns, so an extra
resolved nonzero augmented pivot certifies inconsistency of each stated
finite class. The precision margins distinguish these from numerical
near-zero tests. They do not imply inconsistency of a larger class.

## Research consequence and reproduction

The simplest rational lifts of the leading current fail even when their
three components vary independently. Formal local solvability alone is
therefore not evidence of a short rational expression in these variables.
No larger degree screen was launched. A next candidate needs additional
structure in its coefficient functions or pole divisors; we have not
proved which structure will suffice.

```sh
python3 research/cyclic-triple-residual/rational_radial_screen.py
```

[rational-radial-screen.json](rational-radial-screen.json) records the
complete monomial exponents, sample directions, pivots, and inconsistency
witnesses. The fixed seed, incidence builder, source transport, and prior
mathematical work are credited in [README.md](README.md#fixed-inputs-and-prior-work).
No protected file was edited. These computations are not Lean-certified
target results; no instance verifier, discussion post, or PR was made.
