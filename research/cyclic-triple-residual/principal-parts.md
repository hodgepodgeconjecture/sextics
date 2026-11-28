# A finite critical-locus trace for the principal-parts obstruction

Research checkpoint, 2026-11-28, by OpenAI Codex for
`hodgepodgeconjecture`. This derives an explicit finite-algebra target for
one global obstruction. It does **not** evaluate that trace, finish the
boundary classification, or prove the protected certificate.

## The two obligations in a principal-parts test

Work on the smooth compact normalization of the connected Kummer curve
from [the field audit](incidence-irreducibility.md), over the parameter
field k=C(K,L). Write alpha for the transformed source differential.
If any residue of alpha is nonzero, it cannot be exact. Otherwise each
pole has a local primitive principal part p_P, uniquely determined modulo
regular functions. A rational primitive requires two distinct conditions:

1. A global rational function G realizes all p_P and has no other poles.
2. The resulting holomorphic differential alpha-dG is zero.

The first condition alone is not sufficient. For example, a nonzero
holomorphic differential has empty principal parts and still is not exact.
By the Mittag-Leffler/Serre-duality criterion, the first condition is
equivalent to the vanishing of the sums of residues of p_P theta for every
holomorphic differential theta. This is a standard duality criterion,
not a newly assumed existence theorem for this family.

More explicitly, if
\(D_\alpha=\sum_P\max(0,-\operatorname{ord}_P\alpha-1)P\),
then the exactness question is a finite linear system in the appropriate
character subspace of \(L(D_\alpha)\), imposing \(dG=\alpha\).
The pole divisor still needs to be determined or bounded rigorously.
On a geometrically connected curve, a primitive in the nontrivial
character is unique: the derivative kernel consists of constants, whose
character is trivial.

## One opposite-character pairing

Retain \(r^6=\Delta\), \(w=w_0\), and dots for deformation at fixed
\(\xi,K,L\). Consider

\[
\alpha=\left(\frac{2K^3\xi^3}{3w^3}
             +\frac{\dot\Delta}{18w\Delta}\right)\frac{d\xi}{r},
\qquad \theta=\frac{d\xi}{wr^5}.
\]

The candidate theta is regular at simple discriminant zeros and at the
boundaries previously audited. Its regularity at *all* remaining points
is not yet established. The global obstruction statement below is
conditional on theta being holomorphic.

Put

\[
g=-\frac{\dot\Delta}{3w\Delta_\xi},\qquad
\chi=\frac{g}{r}\theta
=-\frac{\dot\Delta}{3w^2\Delta_\xi\Delta}\,d\xi.
\]

Although defined using the Kummer character, chi is a rational differential
on the degree-120 incidence curve. The exact subtraction is

\[
\alpha-d(g/r)=\left[
\frac{2K^3\xi^3}{3w^3}
+\partial_\xi\left(\frac{\dot\Delta}{3w\Delta_\xi}\right)
\right]\frac{d\xi}{r}.
\]

This is the earlier discriminant-pole subtraction transported to the
endpoint-curve coordinates. Its new use here is to transfer a global
principal-parts pairing to the artificial poles at critical points of
Delta. It is not a new primitive construction.

At an ordinary simple discriminant zero, with w and Delta_xi units,
g/r has the required principal part. The Kummer residue of its product
with theta is

\[
-\frac{2\dot\Delta}{w^2\Delta_\xi^2}.
\]

At an ordinary endpoint branch w=0 with Delta and its w derivative units,
the required principal part is the previously known
\(K^3/(18\xi^2wr)\). The sum over the six unramified Kummer points of its
product with theta is \(-K^3/(36\xi^7\Delta)\). This agrees with
six times the incidence-curve residue of chi: use
\(d\xi=-w\,dw/(12\xi^5)\) and
\(\dot\Delta/\Delta_\xi\to-K^3/(6\xi^2)\).

At a simple critical point of Delta, away from the endpoint, source,
and discriminant singularities, alpha is regular but chi has residue

\[
\operatorname{res}\chi
=-\frac{\dot\Delta}{3w^2\Delta\Delta_{\xi\xi}}.
\]

The global residue theorem therefore transfers the ordinary source-pole
pairing to **positive** contributions
\(2\dot\Delta/(w^2\Delta\Delta_{\xi\xi})\) at these critical points,
plus correction terms at any other, unclassified poles. All sums here
are over geometric points with their residue-field traces, not just the
rational points of a finite-field model.

## Previously audited boundaries contribute zero

For each known boundary family, chi is regular. The existing valuation
bounds give the following orders on the incidence normalization:

| Family | Lower bound for order of chi |
| --- | ---: |
| Infinity, 72 excess-family branches | 15 |
| Infinity, c=0 isolated type | 36 |
| Infinity, b=0 cusp type | 28 |
| Finite constant-cubic family | 24 |
| Finite translated-cubic family | 27 |

For infinity, let ord(xi)=-e, ord(Delta)=m, and let ell bound the pole
order of dot(Delta)/Delta. The bound is \(4e-\ell-m-1\).
At a finite point with ord(xi-xi0)=e and w a unit, it is
\(2e-\ell-m-1\). In the translated-cubic family ell=4, since
Delta_old=h^-6 Delta_raw and h_s/h can have order -4.
These bounds use the proved orders from [infinity-boundary.md](infinity-boundary.md)
and [finite-boundary.md](finite-boundary.md). They do not cover unnamed
boundary families or higher critical points.

## The trace can use polynomial cofactors

Use the notation of [the cofactor construction](discriminant-chart.md):
seven equations Q=coeff(7H), dependent variables y=(xi,a,b,c,d,e,f),
\(T=\det Q_y\), \(L_3=T\partial_{k_3}\),
\(N_3=L_3(\Delta)\), and \(N_s=T\Delta_s\), where this last source
variation holds k3 fixed. Along the fixed-parameter curve,

\[
\partial_\xi=-\frac{6\xi^5}{wT}L_3,\qquad
\Delta_\xi=-\frac{6\xi^5N_3}{wT}.
\]

On N3=0, the change of source coordinate gives dot(Delta)=Ns/T, and

\[
\Delta_{\xi\xi}=
\frac{36\xi^{10}}{w^2T^2}L_3N_3.
\]

Consequently the critical-point contribution simplifies to

\[
\boxed{\frac{T N_s}{18\xi^{10}\Delta\,L_3N_3}}.
\]

No second implicit source derivative is required. Let Acrit be the
finite algebra over C(K,L) defined by Q=0, N3=0 in the eight variables
(y,k3), localized away from xi*w*T*Delta*(L3N3). Its points are exactly
the ordinary critical points in this open chart. The last localization
makes their augmented Jacobian invertible, so this algebra is reduced
and finite. Their total contribution is the multiplication trace

\[
\operatorname{Tr}_{A_{\rm crit}/k}
 \left(\frac{T N_s}{18\xi^{10}\Delta\,L_3N_3}\right).
\]

Thus one well-defined computation remains, together with the omitted
strata corrections and the global holomorphy check for theta. A nonzero
corrected trace would exclude the one-component primitive. A zero trace
would not establish exactness: it is only one dual pairing, and even
vanishing of every principal-parts pairing leaves the holomorphic
remainder test described above.

## Verification and status

Run `python3 research/cyclic-triple-residual/principal_parts.py`.
The script checks the subtraction, three local residue formulas, the
cofactor trace simplification, all five boundary-order bounds, and the
previous modulo-13 actual-source residue. The results are in
[principal-parts.json](principal-parts.json). The global critical algebra
and its trace have **not** been computed; no pole-bound completeness or
exactness conclusion is claimed. This supplies a concrete elimination
and trace problem rather than a larger guessed rational-current ansatz.

The standard principal-parts and duality arguments are applied to the
incidence/source machinery credited in [README.md](README.md#fixed-inputs-and-prior-work).
No protected files, discussion posts, or PRs are changed.
