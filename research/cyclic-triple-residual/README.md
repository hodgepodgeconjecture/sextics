# Cyclic residual: actual-source obstruction and the next current equation

Research checkpoint, 2026-11-26. Prepared for the contributor by OpenAI Codex.
**Research branch checkpoint, saved to the contributor's fork with permission.**
Discussion posts and PRs require the contributor's review.
No target certificate, Lean verification, unrestricted nonexistence theorem,
or acceptance is claimed.

## Result and scope

On the smooth local component through the accepted weighted-line seed, the
actual `cyclic-triple-residual` source does **not** admit a primitive using
only differentiation along the lifted curve, with the endpoint parameters
held constant. This is an all-degree obstruction for that precise normal
representative and family, not a cutoff search.

The first nonzero degeneration term has a nonzero residue at infinity and
a nonzero elliptic remainder for the explicit direction tested. A universal
five-by-five determinant identity reduces this class of cubic-cover
exactness problems to two scalar obstructions. Parameter derivatives can
cancel both leading terms formally; a rational current for the full family
has not been constructed.

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

For the actual family, a conservative degree bound is `deg J <= 6`,
`deg A <= 6`, and `deg R <= 15`, so `deg s <= 13` suffices for the
complete test. This bound follows from the family, not a guessed cutoff.

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

## Next research task

Apply the two-remainder reduction to the full polynomial `R` over the
incidence field, not just its first jet. Use the connection and exact
corrections above to search for a rational solution of the two-component
current equation; verify reconstruction before
attempting Lean transport and clearing. The family is isotrivial after a
finite cover, but moving marked points and source corrections remain.

Do not repeat a search for a time-only primitive, assume the finite trace
removes a nonzero differential class, or infer global rationality from the
leading beta-antiderivatives. No AWS computation is justified yet.

## Reproduction and validation

Requires Python 3 and SymPy (tested with Python 3.12 and SymPy 1.14.0).

```sh
python3 verify_initial_obstruction.py
python3 verify_cubic_connection.py
python3 line_exactness.py --prime 101 --precision 8 --output result.json
python3 line_exactness.py --prime 109 --precision 10 --output result-109.json
```

The first script verifies the characteristic-zero jet, source derivative,
determinant, complete normal form, exact obstruction, and leading parameter
current with symbolic arithmetic. It is the primary argument.
The connection script verifies eight differential identities, six
flatness identities, and a scaling control. The modular script
independently builds all seven incidence equations,
Hensel-lifts them at the exact coupling, computes first and second implicit
derivatives and all determinant-variation terms, then tests the complete
16-by-14 polynomial system. It also reconstructs the full source from its
two-remainder normal form: the normalized remainder pair is
`(38008933,69168482)` modulo `101^4`, and
`(399742904841,92335547497)` modulo `109^6`. These are adic sample values,
not rational reconstructions of the generic coefficient functions.
Both primes give ranks 14 and 15 after
normalizing `v_p(Z)=1` and `v_p(R)=4`; they are implementation cross-checks,
not substitutes for the characteristic-zero argument.

The seed equation is normalized as `7H=7epsilon`. The source uses the
original numerator, not the accepted companion numerator. The modular
program's report records the actual incidence points for replay.

No Lean toolchain checks or protected verifier runs were performed for
this checkpoint. The repository now requires Lean 4.36.0. No repository
protected files were edited. This is a research-branch checkpoint, not a
solution submission; no discussion post or PR has been made.
