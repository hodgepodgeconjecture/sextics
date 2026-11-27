# Discriminant basis, source connection, and local residues

Research checkpoint, 2026-11-27, prepared by OpenAI Codex for
`hodgepodgeconjecture`. This builds exact finite-algebra data for the
integration problem and identifies local residue criteria. It does not
compute the global differential reduction or prove a target certificate.

## The discriminant generates the generic incidence algebra

Starting from [the three slope equations](algebraic-curve-model.md), form
their affine quotient over the endpoint field
\(K(\xi,w_0)\), where \(K=\mathbb C(k_0,k_2)\). The previous degree
calculation gives total generic length 120. We construct multiplication
matrices without finding or approximating the individual line slopes.

The degree-at-most-12 relation matrix has rank 335 among 455 columns;
the degree-at-most-13 matrix has rank 440 among 560 columns. Both leave
a quotient of dimension 120. The degree-12 basis maps invertibly into
the degree-13 quotient. This supplies multiplication by \(a,b,c\).
The program checks that the matrices commute, annihilate all three
incidence relations, and send the unit to the correct basis monomials.

Next reconstruct \(Z\) and form its discriminant \(\Delta\) in this
algebra. At each of the primes 1009 and 1013, the vectors

\[
1,\Delta,\Delta^2,\ldots,\Delta^{119}
\]

are a basis. The resulting monic degree-120 relation is squarefree and
has nonzero constant term. The report records its exact coefficients
and polynomial expressions for all three slopes in this discriminant
basis. The base endpoints are recorded too; they are different exact
finite-field specializations of the fixed target construction.

These invertible-matrix and squarefreeness certificates hold on a generic
characteristic-zero open set: as in the earlier degree calculation, the
endpoint points lift algebraically with nonzero \(w_0\), and the nonzero
minors remain nonzero after lifting. Thus the discriminant is a primitive
element of the **total generic finite incidence algebra** over the
endpoint field. This does not prove that algebra irreducible or identify
the factor corresponding to our seed component.

In particular there is a generic representation

\[
\mathcal A\simeq K(\xi,w_0)[T]/(P(T)),\qquad \Delta\leftrightarrow T,
\qquad \deg_T P=120.
\]

We have computed the polynomial coefficients at exact specializations,
not expanded the global parameter-dependent coefficients of \(P\).
This statement also does not change the earlier warning that the degree
over \(K(\Delta)\) is a different question.

## The actual source in that basis

The code differentiates the reduced incidence equations in two directions:
the endpoint variable \(\xi\), and the auxiliary deformation \(F-sA\)
at fixed \(\xi,k_0,k_2\). It solves the three slope-variation equations
inside the full quotient algebra. Concretely this is a 360-by-360 exact
linear system at each specialization. Its solutions satisfy all the
linearized equations.

Differentiating the reconstructed cubic yields exact vectors for
\(\Delta_\xi\) and \(\dot\Delta\). With \(r^6=\Delta\), the
transported source is

\[
\widetilde\alpha=V\,d\xi/r,\qquad
V=\frac{2k_0^3\xi^3}{3w_0^3}
  +\frac{\dot\Delta}{18w_0\Delta}.
\]

The report gives \(V\), \(\Delta_\xi\), and \(\dot\Delta\) as
120 coefficients in the discriminant basis at each test fiber. The
source-coordinate-change identity is rechecked with these actual vectors.

For the basis \(e_j=\Delta^j/r\), differentiation is

\[
\partial_\xi e_j=(j-1/6)\Delta^{j-1}\Delta_\xi/r,
\qquad 0\le j<120.
\]

Reduction modulo \(P\) gives the connection matrix \(C\); the report
contains enough data to reconstruct it, including \(\Delta^{-1}\).
An additional product-rule check verifies

\[
M_\Delta'+C M_\Delta-M_\Delta C=M_{\Delta_\xi},
\]

where the derivative of the companion matrix is obtained by
differentiating its defining polynomial relation. This checks reduction
of the top power as well as the lower powers.

Consequently a primitive in this character can be sought through
\(g'+Cg=V\), with coefficients on the endpoint curve. The code evaluates
this system's data exactly at two ordinary fibers. It does **not** solve
the differential equation, reconstruct its global coefficient functions,
or compute its residues at singular fibers. We currently retain all 120
sheets: a solution on the whole algebra would suffice, but failure on
the whole algebra would not automatically exclude the selected component.

## Local residue criteria

There are useful deductions before enumerating every ramification point.
On a normalized incidence-curve branch, let \(u\) be a uniformizer and
\(\operatorname{ord}_u\Delta=m\). Put \(g=\gcd(6,m)\) and
\(e=6/g\). After the Kummer extension, a term \(u^n\,du/r\) has
uniformizer exponent

\[
e(n+1)-m/g-1.
\]

A residue requires this exponent to be \(-1\), or \(6(n+1)=m\).
Therefore **if \(6\nmid m\), every differential in this character has
zero local residue**, regardless of the integer pole order of its
coefficient. Unit factors do not affect the congruence. This applies to
zeros and poles of \(\Delta\). If \(6\mid m\), the character alone
does not force vanishing; the case \(m=0\) is included here.

Coordinate ramification also needs care. Write the original differential
as \(\alpha=(\partial_s f)|_x\,dx\), with
\(x=k_3\) and \(f=J/(3r)\). On a smooth local family chart
\(x=x(u,s)\), the exact chain-rule identity is

\[
\alpha=(f_s x_u-f_u x_s)\,du
 =\partial_s(fx_u)\,du-d(fx_s).
\]

If \(f\) and \(x\) are regular in this family chart, the first
expression is regular even when \(x_u=0\). Thus a pole in an implicit
derivative caused by a ramified coordinate projection need not give a
pole or residue of the differential. This is a conditional local
statement; we have not proved that every point of every degeneration
admits such a regular family chart.

For our two endpoint descriptions this also gives concrete open-set tests:

- At \(w_0=0\), with \(\xi,\Delta,T\) units, the original source
  \(B\,dk_3/r\) is regular. The apparent pole of the transported source
  differs from it by the already displayed exact differential.
- At \(\xi=0\), with \(w_0,\Delta,\det S\) units, the transported
  source is regular; its deformation formulas contain no negative power
  of \(\xi\). The old endpoint pole is exact on this open set.

The remaining residue analysis must inspect the actual divisors with
\(6\mid\operatorname{ord}\Delta\), boundary points, and any family
singularities where these regularity hypotheses fail. Even if every
residue vanishes, compact de Rham cohomology can still obstruct a rational
primitive. No global residue-vanishing or exactness theorem is asserted.

## Reproduction

```sh
python3 -m pip install -r research/cyclic-triple-residual/requirements.txt
python3 research/cyclic-triple-residual/finite_curve_algebra.py
```

All finite-algebra, source-variation, connection product-rule, and local
symbolic checks pass. Computations use exact modular arithmetic; NumPy
integer bounds are small enough to avoid overflow at the stated primes.
[finite-curve-algebra.json](finite-curve-algebra.json) records the exact
coefficients and checks. No floating-point roots or fitted current are
used. Prior mathematical work and the fixed-input provenance are credited
in [README.md](README.md#fixed-inputs-and-prior-work).

No protected files were edited. No target verifier, discussion post,
or PR was run or submitted for this research checkpoint.
