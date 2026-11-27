# The generic integration curve is connected

Research checkpoint, 2026-11-27, prepared by OpenAI Codex for
`hodgepodgeconjecture`. This resolves the component ambiguity left by
[the finite-algebra calculation](finite-curve-algebra.md). It does not
prove existence or nonexistence of a rational primitive. The arguments
below are algebraic and reproducibly checked, but not formalized in Lean.

## Irreducibility of the degree-120 incidence field

Let \(B\) be the endpoint threefold over \(\mathbb Q(i)\), with coordinates
\(K,L,\xi,w_0\) and equation

\[
w_0^2=-4\bigl(\xi^6+1+L^6+(1-\kappa^2/4)K^6\bigr),\qquad\kappa=5/7.
\]

Its generic incidence algebra has dimension 120. Multiplication by the
cubic discriminant has a monic degree-120 characteristic polynomial
\(P\) over \(\mathbb Q(i)(B)\), and the discriminant is a cyclic
generator on a nonempty open set. We recompute this algebra at three
good finite-field endpoints. The resulting exact factor constraints are:

| Prime | Irreducible factor degrees, or test | Remaining proper factor degrees of generic \(P\) |
| ---: | --- | --- |
| 1009 | 1, 16, 103 | 1, 16, 17, 103, 104, 119 |
| 1013 | 1, 4, 54, 61 | 1, 119 |
| 1069 | \(\gcd(P,T^{1069}-T)=1\) | none |

The first two rows intersect subset-sum constraints. In the last row,
a degree-119 proper factor would leave a degree-one complementary factor,
so the absence of a linear factor excludes both remaining possibilities.
All three primes split in \(\mathbb Q(i)\).

This argument requires good specialization, not just a factorization at
an arbitrary point. Here the incidence is a relative zero-dimensional
complete intersection of degrees 4, 5, 6, with no points at direction
infinity near each test endpoint. It is finite flat of degree 120 there.
The reconstructed multiplication matrices and cyclic discriminant basis
identify its characteristic polynomial and verify it is squarefree.
The endpoint local rings at these primes are regular (in particular,
\(w_0\ne0\)), hence integrally closed.

If generic \(P\) had a monic proper factor over \(\mathbb Q(i)(B)\),
its coefficients would be integral over each such local ring: they are
elementary symmetric functions of roots of a monic integral polynomial.
They already belong to its fraction field, so normality makes them
regular. The factor therefore specializes monically with the same degree
at every test endpoint. The table excludes that degree. This also
explains why an unspecified denominator in a hypothetical factor cannot
evade the modular tests.

Thus \(P\) is irreducible and the generic incidence algebra is a field.
The script verifies the supplied factorization products and checks each
factor independently by Rabin's finite-field irreducibility criterion;
it does not rely only on reported factor degrees.

## Why this is the selected geometric component

The accepted seed is defined over \(\mathbb Q(i)\):

\[
(\xi,a,b,c,d,e,f;K,L,k_3)=(i,0,0,i,0,0,0;0,0,0).
\]

The seven incidence equations have determinant
\(-135000i/2401\) in the seven dependent variables. The seed is smooth,
with free local parameters \(K,L,k_3\). Since
\(w_0=2k_3+\kappa K^3\) has derivative 2 in \(k_3\), its component
dominates the endpoint threefold. Its first cubic jet is
\(Z=\lambda(\gamma-\alpha t^3)+O(\lambda^2)\), so its discriminant
is not identically zero either. These facts are rechecked from the
protected target.

The irreducible generic incidence field consequently belongs to the
component containing this seed. That component is geometrically integral:
over an algebraic closure, its geometric components are permuted
transitively, while a smooth rational point lies on exactly one of them.
Every conjugation fixes the seed and hence that component, forcing there
to be only one. In particular the incidence degree over
\(\mathbb C(K,L)(\xi,w_0)\) is **exactly 120**, not merely an upper
bound on a possibly smaller selected factor.

This does not assert that the full projective incidence scheme has no
extra components over special endpoint loci.

## The generic curve has no extra algebraic constants

There is a further distinction between geometric integrality of this
threefold and geometric integrality of its generic curve over
\(k=\mathbb C(K,L)\). The latter also holds here.

The constant-cubic family in [the finite-boundary audit](finite-boundary.md)
gives a smooth point of the normalized generic curve with residue field

\[
k\bigl(\sqrt{(\kappa^2-4)K^6-4}\bigr),
\]

using \(\rho=i\). Its field of algebraic constants embeds in this
quadratic residue field: a nonzero element algebraic over \(k\) has
valuation zero at any place trivial on \(k\), and reduction is an
injective field map. Thus a nontrivial constant extension could only be
this particular quadratic extension.

Suppose its square root were in the incidence field after extending
constants to \(\mathbb C\). It is unique up to sign and regular at the
smooth rational seed, where its square is \(-4\) and its value is
\(2i\) or \(-2i\). Every conjugate fixes this value, so the sign is
fixed and the root descends to \(\mathbb Q(i)(M)\).

But at the good fiber modulo 1009 there is a rational incidence point
(the linear factor in the table). At \(K=1\), the proposed square is

\[
-367/49\equiv775\pmod{1009},\qquad 775^{504}\equiv-1.
\]

It is a nonsquare. A rational square root in the function field would be
integral at this point, hence regular because the incidence is etale over
the regular endpoint local ring there. Its reduction would give the
impossible square root in \(\mathbb F_{1009}\). This excludes the
quadratic constant extension. In characteristic zero, having no extra
algebraic constants makes the generic curve geometrically integral.

This assertion concerns the generic curve. It does not claim that every
specialization of \(K,L\) is irreducible or smooth.

## Genuine degree-six Kummer branching

The sixth-root extension is not redundant. An exact witness modulo 13
is the following incidence point, in the chart \(X_1=1,Y_1=t\):

\[
(\xi,a,b,c,d,e,f;K,L,k_3)=(9,10,12,3,6,8,3;1,2,1).
\]

The cubic is \(Z=1+6t+8t^2+3t^3\). Its discriminant is zero,
and its gcd with its derivative has degree one. The seven-equation
Jacobian for \(7H\) has determinant 12 modulo 13. The three implicit
endpoint derivatives of the discriminant are

\[
(\Delta_K,\Delta_L,\Delta_{k_3})=(9,3,6)\pmod{13}.
\]

In particular, the augmented eight-equation Jacobian for incidence plus
\(\Delta=0\), in the dependent variables together with \(k_3\), has
determinant 7. The Jacobian criterion and Hensel lifting give a
characteristic-zero smooth discriminant divisor etale over the
\((K,L)\) parameters near this witness. It is therefore horizontal,
and the discriminant has order one at its generic point on the generic
curve. This is an exact smooth-lifting argument; the finite-field point
alone is not being identified with a complex point.

At that valuation, \(R^6-\Delta\) is Eisenstein. The simple horizontal
zero persists after extending the parameter constants algebraically.
Consequently the Kummer curve is geometrically integral and has degree
6 over the incidence curve, hence **degree 720 over the endpoint curve**.

The actual source is nontrivial at this divisor. The same witness gives
\(w_0=12,\Delta_\xi=4,\dot\Delta=4\) modulo 13, with the dot denoting
the fixed-\(\xi,K,L\) deformation \(H-sA\). Thus the transported source
has a genuine second-kind pole at a generic point of this divisor. A
primitive, if it exists, must have a simple pole there.

For orientation, the local opposite-character pairing contribution is

\[
-\frac{2\dot\Delta}{w_0^2\Delta_\xi^2}\equiv6\pmod{13}.
\]

This is one local contribution, **not the global sum** and not an
obstruction to exactness by itself. Its derivation uses the principal
part \(-\dot\Delta/(3w_0\Delta_\xi r)\) of a possible primitive
and the locally regular differential \(d\xi/(w_0r^5)\).

## Reproduction, credit, and limits

```sh
python3 research/cyclic-triple-residual/finite_boundary.py
python3 research/cyclic-triple-residual/incidence_irreducibility.py
python3 research/cyclic-triple-residual/kummer_branch.py
```

[The irreducibility report](incidence-irreducibility.json) contains the
three fiber polynomials, factor certificates, seed checks, and quadratic
constant-field test. [The branch report](kummer-branch.json) contains the
small discriminant witness and all derivative checks. The incidence and
source constructions credited in [README.md](README.md#fixed-inputs-and-prior-work)
remain the starting point; the seed is due to `rbajaj5`.

These results make the integration field precise and remove component
ambiguities. They do not establish a global pole classification, compute
the compact de Rham obstruction, construct a rational current, or prove
the protected `HasCertificate` statement. No discussion or PR was posted.

The next decisive calculation is a global principal-parts test on this
connected curve. At an ordinary simple discriminant zero the required
principal part is the one displayed above; at ordinary endpoint branch
points the known exact subtraction supplies the principal part. Infinity
and the two audited finite families allow no poles. The remaining finite
boundary and multiple-zero strata must be resolved before these conditions
give a complete Riemann--Roch problem. Solving that finite problem, or
evaluating a nonzero global residue-duality pairing, would settle this
one-component route. It would not settle the less restrictive equation
allowing all three parameter-current components.
