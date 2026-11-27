# Cyclic residual: actual-source obstruction and the next current equation

Research checkpoint, 2026-11-26. Prepared for the contributor by OpenAI Codex.
**Research branch checkpoint, saved to the contributor's fork with permission.**
Discussion posts and PRs require the contributor's review.
No target certificate, unrestricted nonexistence theorem, or acceptance is
claimed. Nine universal reduction lemmas now compile and pass independent
kernel replay under the pinned Lean 4.36.0 toolchain; target binding and the
scalar primitive remain open.

**Strategy reassessment:** [the new geometry audit](strategy-reassessment.md)
explains why connection trivialization is not an existence argument for
the scalar current. It checks the singular root-cover model at the actual
seed, computes the conditional canonical exponents of the line and root
spaces, and screens two alternative quotients by their characters.
Unstructured enlargement of the bounded current ansatz is paused pending
a geometric mechanism or a justified pole bound. Reproduce the six new
exact checks with `python3 verify_strategy_geometry.py` from this folder.

**Paired-cubic follow-up:** [the exchange audit](paired-cubics.md) shows
that the natural quotient exchange annihilates the target character.
The two discriminants have exceptional-divisor orders 4 and 1 at the
seed, excluding a rational horizontal identification of their compact
connections. Five exact checks are reproducible with
`python3 verify_paired_cubics.py`. These are obstructions to the specified
shortcuts, not to the target certificate or arbitrary parameter currents.

**Constructive continuation:** [parameter-current construction](constructive-current.md)
gives an exact rational subtraction removing the explicit discriminant pole
on a suitable chart, and a formal radial solution near the seed. Four exact
finite-field runs verify the latter modulo the twentieth power of the
radial parameter; its leading current is checked symbolically in
characteristic zero. Rationality of a full current remains unproved.

## Result and scope

On the smooth local component through the accepted weighted-line seed, the
actual `cyclic-triple-residual` source does **not** admit a primitive using
only differentiation along the lifted curve, with the endpoint parameters
held constant. This is an all-degree obstruction for that precise normal
representative and family, not a cutoff search.

The first nonzero degeneration term has a nonzero residue at infinity and
a nonzero elliptic remainder for the explicit direction tested. A universal
five-by-five determinant identity reduces this class of cubic-cover
exactness problems to two scalar obstructions.

**Latest progress:** an exact change of source representative now removes
the entire residue component, using three explicit rational parameter
currents on the full incidence field. Its time numerator has degree three
instead of nine. The remaining equation is scalar, with connection
`d - (1/6)dlog(Delta)`. Its compact elliptic component is provably nonzero;
a rational parameter current for that remaining component has not been
constructed. The preceding time-only obstruction remains valid for its
stated representative. No auxiliary deformation changes the fixed target.

The known quotient certificate does not extend to the residual target just
by allowing rational functions in its integration step.

## Fixed inputs and prior work

Upstream commit: `d22d25556ed2de562945e3fb78d1b1cec489f367`.
The unchanged record is included under `inputs/`; its source is
`A = X0^3 Y0^3`, with `q = 1` and coupling exactly `5/7`.

Primary sources:

- [Exact target](https://github.com/hodge-computer/sextics/blob/d22d25556ed2de562945e3fb78d1b1cec489f367/data/cyclic-triple-residual/instance.rec).
- [Accepted companion construction](https://github.com/hodge-computer/sextics/blob/d22d25556ed2de562945e3fb78d1b1cec489f367/data/cyclic-triple-fano/solution.md).
- [Residual discussion, read in the preceding audit](https://github.com/hodge-computer/sextics/discussions/155).
- [Curve exactness criterion, book section 5.5](https://github.com/hodge-computer/sextics/blob/d22d25556ed2de562945e3fb78d1b1cec489f367/docs/05-computation.md).

The seed and its freeness calculation are due to `rbajaj5`. The accepted
incidence construction and source machinery credit `fqlx`, `HeadOfNaps`,
Benjamin E. Diamond, and the Lean/mathlib contributors. The prior residual
thread contains `fqlx`'s marked-connection/scalar investigations,
`kejace`'s degree-at-most-six exclusions, and `benediamond`'s analysis of
lifted lines, their isotrivial elliptic interpretation, and the excluded
double-root correspondence. This checkpoint uses the anchored family from
the accepted companion and supplies a direct actual-source test on that
family. It makes no novelty claim relative to unpublished or unreplayed
attachment contents.

## Actual-source transport

Use quotient coordinate order `(Y0,Y1,X0,X1,Y2,Z)` and set `Z=X2^3`.
Let `H` be the weighted sextic in the accepted explanation. Enlarge the
anchored incidence by the normal parameter epsilon:

\[
H(\phi(t;\epsilon,k))=\epsilon,\quad
\phi=(\xi+at,t,k_0+bt,k_1,k_2+ct,k_3+dt+et^2+ft^3).
\]

There are seven coefficient equations for `(xi,a,b,c,d,e,f)` and four
independent endpoint parameters `k`. At the accepted seed their Jacobian
is invertible. Implicit differentiation therefore defines the first and
second normal/endpoint derivatives used below. Put

\[
J=\det[\phi_\epsilon,\phi_t,\phi_{k_0},\phi_{k_1},
       \phi_{k_2},\phi_{k_3}],\qquad w^3=Z(t).
\]

The physical volume density acquires the factor `1/(3w^2)`. For the
double pole, its residue source coefficient is the normal derivative

\[
\partial_\epsilon\left(\frac{AJ}{3w^2}\right)
 =\frac{R}{9w^5},\qquad
R=3Z\,\partial_\epsilon(AJ)-2AJ\,\partial_\epsilon Z.
\]

Thus both the derivative of the moving determinant and that of the cover
coordinate are retained. There is a fixed overall orientation sign when
returning to the repository coordinate order; it does not affect any
zero/nonzero or exactness assertion here. A certificate assembly would
have to bind that sign explicitly.

## A complete polynomial test

Assume characteristic zero and that `Z` is a squarefree cubic. The source
differential has deck character `w -> zeta w` acting by `zeta`. Projection
of any rational primitive onto this character gives a primitive of the
form `w r(t)`. At a finite branch point, the differential has pole order
at most three in the uniformizer `w`, so a primitive has pole order at
most two. It has no other finite poles. Consequently

\[
u=\frac{s(t)}{3w^2},\qquad s\in E[t],\qquad
du=\frac{3Zs'-2Z's}{9w^5}\,dt.
\]

The condition is exactly `L_Z(s)=R`, where `L_Z(s)=3Zs'-2Z's`.
No finite algebraic extension of the curve field can repair failure:
normalized trace would descend any primitive.

If `Z=a t^3+b t^2+c t+d`, then on a polynomial of degree `n` the
leading multiplier is `3a(n-2)`. Terms of degree at least five in `R`
can therefore be removed uniquely using `L_Z(t^n)` for `n>=3`.
For the remaining degree-at-most-four problem,

\[
\det[\operatorname{coeff}(L_Z(1)),\operatorname{coeff}(L_Z(t)),
\operatorname{coeff}(L_Z(t^2)),\operatorname{coeff}(Z),
\operatorname{coeff}(tZ)] =8a\operatorname{disc}(Z).
\]

The coefficient rows are degrees zero through four. Hence every polynomial
`R` has a unique decomposition

\[
R=L_Z(s)+C_0 Z+C_1 tZ,
\quad
\frac{R\,dt}{9w^5}
=d\left(\frac{s}{3w^2}\right)
 +\frac{C_0}{9}\frac{dt}{w^2}
 +\frac{C_1}{9}\frac{t\,dt}{w^2}.
\]

The first remaining differential is holomorphic and nonzero; the second
has nonzero residues at the three points at infinity. Exactness holds
if and only if both coefficients vanish. Finite branch points have zero
residues automatically by their local exponents modulo three.

The constant-determinant lemma below sharpens the actual family's bounds
to `deg R <= 9` and `deg s <= 7`. The original checkpoint used the valid
but looser bounds 15 and 13. These bounds follow from the family, not a
guessed cutoff.

## The moving determinant is constant

After expanding along the `Y1/time` and `X1/k1` entries, the evaluation
determinant is the determinant of a four-by-four polynomial matrix `M`
with rows `(Y0,X0,Y2,Z)` and columns `(epsilon,k0,k2,k3)`. Differentiating
`H(phi)=epsilon` gives

\[
(H_{Y_0},H_{X_0},H_{Y_2},H_Z)M=(1,0,0,0).
\]

Multiply by the adjugate. Then `H_Y0 * det(M)` is the three-by-three
minor in the last three rows and columns. Its degree is at most
`1+1+3=5`. Generically `H_Y0(phi)` has degree five, while `det(M)` is
polynomial, so `det(M)` is constant in time. The required degree-five
coefficient is nonzero on our component: along the first jet above its
leading term is `-36 beta^2 lambda^2/kappa`, when `beta != 0`.

Evaluation at `t=0` gives `J=1/(6xi^5)`. This is an identity on the
generic incidence field and hence on its domain of definition, including
specializations where the degree-five coefficient vanishes. It holds
with epsilon variable too. The endpoint equation gives
`xi_epsilon=1/(6xi^5)`, so

\[
J=\frac1{6\xi^5},\qquad J_\epsilon=-\frac5{36\xi^{11}}.
\]

The same argument applies in the `X0`-anchored chart below, using rows
`(Y0,X1,Y2,Z)`. There the degree-five gradient coefficient is nonzero at
the transformed limiting line. This gives a cheaper actual-source
evaluator requiring only the first implicit normal derivative. The code
retains its separate second-derivative determinant computation and checks
agreement, instead of using the simplification as its own verification.

## Characteristic-zero first jet

Keep a symbolic nonzero coupling `kappa`. At the seed with endpoint
`(0,1,0,0)`, allow `Y0=r` and `epsilon=r^6+1`. Perturb the endpoints to

\[
(k_0,k_1,k_2,k_3)=(\alpha\lambda,1,\beta\lambda,\gamma\lambda).
\]

The coefficient Jacobian gives, to first order in lambda,

\[
\begin{aligned}
Y_0&=r-\frac{\kappa\gamma}{6r^4}\lambda t+O(\lambda^2),&
X_0&=\lambda\left(\alpha-\frac{6\beta}{\kappa}t\right)+O(\lambda^2),\\
Y_2&=it+\beta\lambda+O(\lambda^2),&
Z&=\lambda\left(\gamma-\frac{i\alpha}{r}t^3\right)+O(\lambda^2).
\end{aligned}
\]

Here `Y1=t`, `X1=1`. Direct substitution verifies the first-order
incidence equations as identities in `r`, not merely at `r=i`.
At lambda zero, `J=1/(6r^5)` and the normal derivative is
`(6r^5)^(-1) partial_r`. Evaluating at `r=i` gives

\[
Z=\lambda Z_1+O(\lambda^2),\quad Z_1=\gamma-\alpha t^3,
\quad R=\lambda^4 R_4+O(\lambda^5),
\]
\[
R_4=\left(\frac{\alpha t^3}{9}-\frac{\gamma}{6}\right)
     \left(\alpha-\frac{6\beta}{\kappa}t\right)^3.
\]

For `alpha=1,beta=2,gamma=3,kappa=5/7`, the coefficient of `t^4`
in `R4` is `-28/5`. It is zero in `L_(3-t^3)(s)` for every polynomial
`s` of degree at most four; higher-degree polynomials cannot produce a
right side of degree six. Equivalently, the exact rational system has
rank 5 and augmented rank 6. On `w^3=3-t^3`, the residue of
`R4 dt/(9w^5)` at an infinity point with `w/t -> c`, `c^3=-1`, is
`28/(45c^5)`, which is nonzero.

More fully, the two remainder coefficients are

\[
C_0=-\frac{5\alpha^4\kappa^3-216\beta^3\gamma}
              {36\alpha\kappa^3},\qquad
C_1=\frac{2\alpha^2\beta}{\kappa}.
\]

They cannot vanish simultaneously when `alpha*gamma*kappa != 0`:
`C1=0` forces `beta=0`, giving `C0=-5alpha^3/36 != 0`.
In the displayed direction, they are `1777487/4500` and `28/5`.

To justify the passage from a jet to the generic obstruction, use the
uniform finite matrix above. Its leading operator for `Z1` is injective.
A Laurent-series solution would have to start at lambda order three;
its initial coefficient would solve `L_Z1(s3)=R4`, which is impossible.
Alternatively, consistency of the generic finite system would make its
augmented minors vanish identically, including along this formal arc.
This argument does not specialize a primitive at a pole of its unknown
parameter denominator. The invertible incidence Jacobian places the
formal arc on the smooth local component through the accepted seed.

## What parameter currents already do at leading order

At this leading order `Z1` is independent of beta, and

\[
H_0=-\frac{5\alpha^3\beta}{36}
       +\frac{3\beta^4\gamma}{2\alpha\kappa^3},\qquad
H_1=\frac{\alpha^2\beta^2}{\kappa}
\]

satisfy `partial_beta H0=C0` and `partial_beta H1=C1`. Thus a beta
current cancels both leading remainder coefficients in the relative
differential calculation. In coordinates with `alpha=1`,
`k0=lambda,k2=beta*lambda,k3=gamma*lambda`, the endpoint volume change
has a beta-independent factor `lambda^2`, so it does not spoil this
leading cancellation. The sign of a differential-form current depends
on the chosen wedge ordering.

This is **not** a global primitive: higher-order `Z` depends on beta;
its derivative contributes additional terms. Formal power-series
integration also does not establish rationality in the actual incidence
field. These are the obligations for continuation.

## Exact parameter connection and the remaining equation

The calculation now also includes the universal parameter connection,
with exact time corrections. Write `Z=a t^3+b t^2+c t+d`,
`Delta=disc(Z)`, and `omega=(dt/w^2,t dt/w^2)`. Use the column convention
`partial_v omega_j = sum_i Gamma_v[i,j] omega_i + d_t(S_vj/w^2)`.
Then

\[
\Gamma=\begin{pmatrix}
-\tfrac16d\log\Delta&\eta\\
0&-\tfrac23d\log a
\end{pmatrix},
\]
\[
\eta=\frac{1}{6\Delta}\left(
\frac{b\Delta_a+c\Delta_b+d\Delta_c}{a}\,da
-\Delta_a\,db-\Delta_b\,dc-\Delta_c\,dd\right).
\]

The eight rational corrections `S_vj` are supplied in `connection.json`.
They are computed by the same five-by-five inverse, and verified through
the literal polynomial identity

\[
-2t^j\partial_v Z=L_Z(S_{vj})
 +3Z(\Gamma_v[0,j]+t\Gamma_v[1,j]),\qquad j=0,1.
\]

All six curvature identities `partial_u Gamma_v-partial_v Gamma_u +
[Gamma_u,Gamma_v]=0` and the scalar-cubic scaling identity were checked
symbolically. This retains the coupling between the residue direction and
the holomorphic direction; diagonal rescaling alone would omit `eta`.

Let `C=(C0,C1)` now mean the remainders of the **full** source over the
actual finite incidence field `E`, and pull this connection back along
its cubic coefficient functions `a,b,c,d`. A sufficient finishing equation
is to find four rational columns `u_i in E^2` satisfying

\[
\sum_{i=0}^3(\partial_{k_i}u_i+\Gamma_{k_i}u_i)=\frac19 C.
\]

The exact corrections `S_vj` then reconstruct the omitted time component
of the current. This equation is not solved here. It is a sufficient
ansatz with endpoint coefficients and these two relative classes;
additional moving-pole or residue currents remain available outside it.
Finite trace and legal homogeneous clearing would still be required for
the protected certificate.

## Earlier two-component research task

The earlier approach was to apply the two-remainder reduction to the full polynomial `R` over the
incidence field, not just its first jet. Use the connection and exact
corrections above to search for a rational solution of the two-component
current equation; verify reconstruction before
attempting Lean transport and clearing. The source-deformation reduction
below supersedes this with a scalar sufficient equation. The family is
isotrivial after a finite cover; that fact alone does not integrate the source.

Do not repeat a search for a time-only primitive, assume the finite trace
removes a nonzero differential class, or infer global rationality from the
leading beta-antiderivatives. No AWS computation is justified yet.

## Continuation: source-zero infinity and Euler audit

The next test changes the parameterization so the point at infinity has
`X0=0`, where the target numerator vanishes. Let `h=b/k0` for the old
`X0=k0+b t`, and introduce `T=t/(1+h t)`. With `g=1-h T`, scale ordinary
coordinates by `g` and `Z` by `g^3`. The new coefficients are

\[
\begin{aligned}
a'&=a-\xi h,&b'&=-k_1h,&c'&=c-k_2h,\\
d'&=d-3k_3h,&e'&=e-2dh+3k_3h^2,&
f'&=f-eh+dh^2-k_3h^3.
\end{aligned}
\]

Now `X0=k0` and `X1=k1+b'T`; the endpoint is unchanged. All six
coordinate identities and both charts' coefficient equations are checked
against the copied exact `instance.rec`. The fixed polynomial and source
are not changed. We solve the new normal equation `H=epsilon`, rather
than pretending that the old normal simply scales: away from `epsilon=0`
the coordinate transformation sends `H` to `g^6 H`.

The `source_root_chart.py` calculation tracks absolute p-adic precision
through inversions with negative valuation. At endpoints
`(k0,k1,k2,k3)=(p,1,2p,3p)`, **both** remainders are nonzero, with
valuation two at `p=101` and `p=109`. The outputs certify the incidence,
inverse, source determinant, and reconstructed normal form to the
reported precisions (60 input digits at 101, 40 at 109). In this chart
`deg R <= 6`, so the time-primitive bound is only four. Moving infinity
to a zero of the source therefore does not by itself cancel its residue
for this rebuilt normal representative. A claim that the change of chart
finishes the problem would omit parameter and normal-change currents.

There is also an exact scaling obstruction to a naive Euler shortcut.
The endpoint Euler field is
`S=k0*d_k0+k1*d_k1+k2*d_k2+3k3*d_k3`, with divergence six.
The cubic coefficient weights `(a,b,c,d)` in the connection formula are
`(0,1,2,3)`, so `Gamma_S=diag(-1,0)`. Uniqueness of the homogeneous
normal form gives source remainder weights `(-5,-6)`. Therefore

\[
S(C)+(\operatorname{div}S)C+\Gamma_S C=0.
\]

The source is precisely at the zero Euler weight: division by a nonzero
weight cannot supply a primitive. Neither the chart test nor this scaling
identity excludes general rational parameter currents. The continuation
below removes one component by an explicit source-transport identity.

## Exact source deformation removes the residue component

Introduce an auxiliary scalar `s` and solve the same anchored incidence for

\[
H(\phi_s(h,k,t))-s A(\phi_s(h,k,t))=h,
\qquad A=X_0^3Y_0^3.
\]

All expressions below are evaluated at **`s=h=0`**. Thus `F,q,A` remain
exactly the served target. The auxiliary derivative is a way of computing
rational coefficients in the existing incidence field, not a new permitted
certificate differentiation direction. The seven coefficient equations
give `Q_y U=7 coeff_t(A)` for `U=partial_s phi_s`; the already nonsingular
Jacobian makes `U` rational in that field.

Here is an explicit transport current. Omit the fixed `Y1=t,X1=k1` rows
and put `x=(Y0,X0,Y2,Z)`, `z=(h,k0,k2,k3)`,
`M=partial_z x`, `J=det M`, and `V=partial_h x`. For
`i in {0,2,3}`, let `S_i` be the determinant obtained by replacing the
`k_i` column of `M` with `U`. The `k1` current is zero. Since
`dH(U)=A`, `dH(V)=1`, and `dH(partial_ki x)=0`, Cramer's rule gives

\[
JU=AJV+\sum_{i\in\{0,2,3\}}S_i\,\partial_{k_i}x.
\]

The divergence of a row of the cofactor matrix of a Jacobian is zero
(the classical Piola identity). Applying it to this four-variable map gives

\[
J_s=\partial_h(AJ)+\sum_i\partial_{k_i}S_i.
\]

The `Z` row of Cramer's identity also gives
`J Z_s=AJ Z_h+sum_i S_i Z_ki`. Combining these identities, with `w^3=Z`,
proves the **literal rational identity**

\[
\boxed{\quad
\partial_h\!\left(\frac{AJ}{3w^2}\right)
=\partial_s\!\left(\frac{J}{3w^2}\right)
 -\sum_i\partial_{k_i}\!\left(\frac{S_i}{3w^2}\right).
\quad}
\]

Consequently the currents `-S_i/(3w^2)` transport the full original source
to the new one. These currents are rational on the full component;
they are not truncated formal antiderivatives. Each `S_i` is a polynomial
of degree at most six in `t` because the three ordinary rows of `M,U`
are linear and their weighted `Z` row is cubic. The cofactor identity is
verified over characteristic zero with arbitrary Jacobian entries and
arbitrary symmetric Hessians in `verify_source_deformation.py`.

The constant-determinant argument continues to apply to `H-sA`, which is
weighted homogeneous of degree six. Its endpoint equation is

\[
\xi^6+k_0^6+k_1^6+k_2^6+k_3^2+\kappa k_0^3k_3-s k_0^3\xi^3=h.
\]

Thus, at `s=h=0`,

\[
J=\frac1{6\xi^5},\qquad
\xi_s=\frac{k_0^3}{6\xi^2},\qquad
J_s=-\frac{k_0^3}{18\xi^8}.
\]

Both `J` and `J_s` are constant in `t`, while `Z_s` is cubic. The new
source is therefore

\[
\widetilde\omega=\frac{\widetilde R}{9w^5}\,dt,
\qquad \widetilde R=3ZJ_s-2JZ_s,\qquad
\deg_t\widetilde R\leq3.
\]

In the complete five-column normal form, only the `tZ` column has a
degree-four term. Hence its coefficient is **identically zero** for this
new source. Equivalently, it has no residues at infinity. The connection
formula determines the remaining coefficient exactly:

\[
\widetilde\omega=d_t\!\left(\frac{J S_{s0}}{3w^2}\right)
 +B\,\frac{dt}{w^2},\qquad
B=\frac{J_s}{3}-\frac{J\Delta_s}{18\Delta}.
\]

Here `S_s0=sum_v (partial_s v) S_v0`, for the four cubic coefficients
`v`, using the exact corrections in `connection.json`. In the previous
normalization, `C0=9B=3J_s-J Delta_s/(2 Delta)` and `C1=0`.

This is genuine cancellation of the residue component by explicit
parameter currents, rather than a claim that the old residue was zero.
It uses standard determinant transport and Gauss-Manin differentiation;
no novelty claim is made for those general identities.

The same argument proves a broader statement for this incidence chart:
**every weighted-degree-six polynomial source on the quotient has a
representative in the compact elliptic line**, provided one works on the
same open set where the incidence Jacobian, the cubic leading coefficient,
and its discriminant are nonzero. Replace `A` by that polynomial in the
auxiliary deformation. Weighted homogeneity preserves the determinant
degree argument, so `J_s` is still time-constant and the new numerator
still has degree at most three. This generalizes the residue elimination;
it does not prove rational exactness of the compact component for any of
those sources. In particular, it supplies no new target certificate by itself.

### The compact component remains nonzero

An independent characteristic-zero calculation expands the original
incidence to order `lambda^2`, the source deformation to order `lambda^4`,
and the discriminant derivative to order `lambda^7`. For the same endpoints
`k=(alpha*lambda,1,beta*lambda,gamma*lambda)` and branch `xi -> i`, it gives

\[
\begin{aligned}
C_0={}&\lambda^3\left[-\frac{5\alpha^3}{36}
 +\beta^2\left(\frac{24\alpha}{\kappa^4}
               -\frac{6i\alpha^2}{\gamma\kappa^3}\right)\right.\\
&\left.\hspace{24mm}
 +\beta^3\left(\frac{6\alpha}{\gamma\kappa^3}
               +\frac{96i}{\kappa^4}
               +\frac{9\gamma}{\alpha\kappa^3}\right)\right]
 +O(\lambda^4).
\end{aligned}
\]

At `alpha=1,beta=2,gamma=3,kappa=5/7`, the bracket is
`22618411/22500 + (1830248/625)i`, which is nonzero. At `beta=0` it is
`-5alpha^3/36`. Thus the new representative still has a compact elliptic
obstruction to time-only integration. A full certificate has not appeared
merely by removing the residue term.

The independent adic evaluator reconstructs the actual original source,
all three determinant currents and their endpoint derivatives, and the
new source. It checks the full transport identity to at least 30 digits
at both 101 and 109. The exact leading coefficient above reduces to
18 and 65 respectively, agreeing with the independent adic values
(the chosen images of `i` are 10 and 33). This is a cross-check of the
characteristic-zero derivation, not an inference of an identity from samples.

### The remaining scalar equation

It now suffices to find four rational functions `b_i` in the actual
incidence field satisfying

\[
\boxed{\quad
\sum_{i=0}^3\left(\partial_{k_i}b_i
 -\frac{\partial_{k_i}\Delta}{6\Delta}b_i\right)=B.
\quad}
\]

The holomorphic line is preserved by the connection, so these currents
do not reintroduce the residue component. Its exact time corrections,
together with the three already constructed currents, would reconstruct
the original source. This scalar equation is **unsolved** here.

One useful equivalent formulation is obtained on the finite cover
`delta^6=Delta`: the equation becomes
`sum_i partial_ki(b_i/delta)=B/delta`. Also
`B/delta=(1/3) partial_s(J/delta)` at `s=0`. This identifies an ordinary
rational divergence problem on a specified algebraic cover. Trivializing
the connection in this way does not prove that the divergence equation
has a rational solution. Formal integration near the seed does not
establish that either. The Euler resonance also persists.

The next task is to exploit this scalar equation and the incidence
relations to construct an actual rational primitive, retaining all
denominators and their hypotheses. Finite trace, homogeneous clearing,
and target-level Lean verification require that missing construction.

## Certified part of the reduction

`lean/CyclicCubicReduction.lean` proves nine universal lemmas over a field
with the stated characteristic and derivation hypotheses:

- `compact_polynomial_identity` verifies the explicit quadratic correction
  numerator for an arbitrary cubic and an arbitrary first coefficient jet.
- `cubic_derivation` and `discriminant_derivation` bind the coefficient jets
  to differentiation of the actual cubic and discriminant polynomials.
- `cubic_primitive_rule` differentiates `P/(3w^2)` in the cubic cover field,
  retaining the covering-coordinate derivative and the nonzero denominator.
- `compact_rational_identity` verifies the residual scalar coefficient.
- `compact_source_of_quadratic` combines those identities with the derivation.
- `compact_source` instantiates the explicit cubic correction. Its conclusion
  is the exact time derivative plus `(J_s/3-J*Delta_s/(18Delta))/w^2`.
- `kummer_derivative` proves the sixth-root connection trivialization.
- `current_reconstruction` proves that a solution of the scalar current
  equation reconstructs the time and parameter currents. **The scalar
  equation is a premise, not a proved existence claim.**

The module compiles without warnings under Lean 4.36.0 and the repository's
mathlib revision `fc4f59e577feb1a4956671bc04ac18aaaf8cd842`. Its declarations
were independently replayed by `leanchecker`. All audited axioms are among
`propext`, `Classical.choice`, and `Quot.sound`. The measured compile and
replay each take about two seconds, with peak child RSS about 2.2 GiB.
The source hash, exact timings and axiom audit are in `lean-checks.json`;
the two text logs are included.

This certifies the universal cubic algebra and the conditional current
assembly. It does **not** yet certify the seven-equation incidence binding,
the determinant transport on that field, an actual scalar primitive,
trace/clearing, or the protected `HasCertificate` statement. There is no
replacement target definition or assumed certificate in the module.

To replay after installing the pinned toolchain, from the repository root:

```sh
cd formal/lean
lake exe cache get Mathlib.RingTheory.Derivation.Basic Mathlib.Tactic.LinearCombination Mathlib.Tactic.FieldSimp Mathlib.Tactic.Ring
cd ../..
python3 research/cyclic-triple-residual/run_lean_checks.py
```

The script also accepts `--formal-root` and `--lake` for an explicit local
installation. It reads the repository's existing options and manifest,
builds temporary artifacts, rejects compilation warnings and disallowed
axioms, and separately invokes the kernel checker. It does not run or
replace the protected instance verifier.

## Bounded scalar-current search: a combined class is excluded

The next search tests the following **specific** finite class, combining
endpoint rational expressions and linear dependence on the weight-zero
curve coefficients:

\[
b_i=J\left(\frac{P_i(\xi,k_0,k_1,k_2,k_3)}{\xi^2}
 +\sum_{\ell\in\{a,b,c,f\}}\ell\,
 Q_{i\ell}(\xi,k_0,k_1,k_2,k_3)\right).
\]

Give `(xi,k0,k1,k2,k3)` weights `(1,1,1,1,3)`. Each `P_i` has
weighted degree `weight(k_i)+2`; each `Q_iell` has degree `weight(k_i)`.
All their coefficients are constants, not functions of the parameters.
There are **261 coefficients**. This class contains both the endpoint-only
denominator-`xi^2` ansatz and the polynomial ansatz linear in `a,b,c,f`,
so the test permits cancellation between these two parts.

The exact Euler zero currents `b_i=J*f*weight(k_i)*k_i` remove a known
14-dimensional kernel: ten choices `f=M_2/xi^2` with `M_2` a weight-two
endpoint monomial, and `f=a,b,c,f`. The Euler field has divergence six,
`S(J)=-5J`, `S(Delta)=6Delta`, and `S(f)=0`, so these are identically
zero under the scalar current operator. Setting the corresponding `k0`
coefficients to zero is an exact choice of representative, leaving
**247 unknowns**; it imposes no restriction on existence within this class.

At each sampled endpoint, the script Hensel-lifts the original seven
incidence equations and computes all four endpoint derivatives and the
source-deformation jet with their actual Jacobian inverse. For a basis
term `b_i=J*M_i/xi^m`, the cleared scalar equation has column

\[
18\Delta\xi\,\partial_{k_i}M_i
-18(m+5)\Delta M_i\,\partial_{k_i}\xi
-3\xi M_i\,\partial_{k_i}\Delta
\]

and right side `-2 Delta k0^3 xi^(m-2)-xi^(m+1) Delta_s`.
These formulas are independently checked by symbolic differentiation.
In the combined class `m=2`, with the polynomial slope terms represented
by numerators divisible by `xi^2`. At every point the code also checks
the Euler derivatives of `xi`, `Delta`, and the four weight-zero slopes.

At both primes, 248 sampled directions yield coefficient rank **247**
and augmented rank **248**. The final transformed right side has
valuation nine, is known modulo `p^44`, and has leading residue 45 at
101 and 56 at 109. Input precision is 48 digits. Precision is propagated
through every elimination and inversion, including zero representatives;
no floating-point rank tolerance or rational reconstruction is used.

This is an inconsistency certificate for this finite coefficient class.
The nonsingular incidence Jacobian supplies actual `Q_p` points behind
the congruences. The 247 pivots and the resolved nonzero augmented pivot
prove the rank claims over `Q_p`. An identity on the seed component with
constant coefficients in `Q(i)` would evaluate at all these points and
is therefore impossible. The exclusion also descends from complex
constant coefficients: expansion at the seed gives a linear system over
`Q(i)` in finitely many unknown coefficients, and any consistent complex
solution of that system would admit a solution over `Q(i)`.

The scripts and results are **not Lean-verified**. In particular, this
does not exclude currents with nonlinear slope dependence, dependence on
the other incidence coefficients, larger endpoint denominators, poles on
the incidence Jacobian or discriminant, or other supports. It is not an
unrestricted nonexistence result. No further blind increase of these
bounds is justified by this calculation alone.

Reproduce the combined screen with:

```sh
python3 scalar_current_screen.py --prime 101 --precision 48 --mixed --sample-limit 270 --output scalar-mixed-101.json
python3 scalar_current_screen.py --prime 109 --precision 48 --mixed --sample-limit 270 --output scalar-mixed-109.json
```

Each run takes about 90 seconds on the current workspace and stops at the
first resolved inconsistent row. The compact reports give the random seed,
sample counts, pivot bounds, and hashes of the reproducible directions and
pivot records. Add `--full-report` to retain those detailed records. The next
mathematical step must use more of the incidence field or find a serving
correspondence; the scalar equation and target proof remain unresolved.

## Reproduction and validation

Requires Python 3 and SymPy (tested with Python 3.12 and SymPy 1.14.0).

```sh
python3 run_checks.py
# Individual stages:
python3 verify_initial_obstruction.py
python3 verify_cubic_connection.py
python3 verify_chart_and_determinant.py
python3 line_exactness.py --prime 101 --precision 8 --output result.json
python3 line_exactness.py --prime 109 --precision 10 --output result-109.json
python3 source_root_chart.py --prime 101 --precision 60 --output source-root-101.json
python3 source_root_chart.py --prime 109 --precision 40 --output source-root-109.json
python3 verify_source_deformation.py
python3 source_deformation.py --prime 101 --precision 30 --output source-deformation-101.json
python3 source_deformation.py --prime 109 --precision 30 --output source-deformation-109.json
```

The first script verifies the characteristic-zero jet, source derivative,
determinant, complete normal form, exact obstruction, and leading parameter
current with symbolic arithmetic. It is the primary argument.
The connection script verifies eight differential identities, six
flatness identities, and a scaling control. The modular script
independently builds all seven incidence equations,
Hensel-lifts them at the exact coupling, computes first and second implicit
derivatives and all determinant-variation terms, then tests the complete
10-by-8 polynomial system (tightened from the initial 16-by-14 system).
It also reconstructs the full source from its
two-remainder normal form: the normalized remainder pair is
`(38008933,69168482)` modulo `101^4`, and
`(399742904841,92335547497)` modulo `109^6`. These are adic sample values,
not rational reconstructions of the generic coefficient functions.
Both primes give ranks 8 and 9 after
normalizing `v_p(Z)=1` and `v_p(R)=4`; they are implementation cross-checks,
not substitutes for the characteristic-zero argument.

The seed equation is normalized as `7H=7epsilon`. The source uses the
original numerator, not the accepted companion numerator. The modular
program's report records the actual incidence points for replay.

The Python checks are distinct from the Lean checks above. No protected
instance verifier run has been made, since there is no complete target
proof. No repository protected files were edited. This is a research-branch
checkpoint, not a solution submission; no discussion post or PR has been made.
