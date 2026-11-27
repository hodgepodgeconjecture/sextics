# Reassessing the scalar-current route

Research checkpoint, 2026-11-26, prepared by OpenAI Codex for
`hodgepodgeconjecture`. This is a research note, not a solution submission.
The contributor has authorized saving to this branch; discussion posts
and PRs await their review.

## Decision

Pause enlargement of the bounded rational-current ansatz. Keep the exact
source transport and the nine checked universal Lean lemmas, but do not
interpret their existence as evidence that the missing rational primitive
is close. They simplify a necessary piece of the construction without
establishing an existence mechanism for that piece.

The next useful work is to identify an algebraic correspondence that serves
the target character, with a checkable nonzero action on this particular
source. A new candidate must clear that character test before we invest in
its incidence algebra or formalization. Nothing below proves that the
current route cannot work, or that the target has no certificate.

No AWS computation is justified by the present bottleneck.

## What the reduction actually achieved

With the notation of [the research account](README.md), the remaining
sufficient equation on the incidence field is

\[
\sum_i\left(\partial_{k_i}b_i-
 \frac{\partial_{k_i}\Delta}{6\Delta}b_i\right)
 =B,\qquad
B=\frac{J_s}{3}-\frac{J\Delta_s}{18\Delta}.
\]

On adjoining \(\delta^6=\Delta\), this becomes

\[
\sum_i\partial_{k_i}(b_i/\delta)
 =B/\delta=\tfrac13\partial_s(J/\delta).
\]

This is a change of presentation of a differential equation on an
algebraic field. It is not integration on a rational affine space. There
is no theorem here saying that the right side is a rational divergence.
Even on a rational field, arbitrary rational differential forms need not
be exact. Finite extensions do not automatically make them exact either.

The previous modular rank failures exclude the explicitly stated
261-coefficient class only. They supply neither an unrestricted
nonexistence theorem nor a principled next degree bound. The Lean
`current_reconstruction` lemma assumes the scalar current as a premise;
that premise remains the main missing mathematics. The exact fixed
`HasCertificate` statement has not been proved or submitted to its
protected verifier.

## The line space is not supplied with rational geometry

Write the quotient equation, with the fixed coupling \(\kappa=5/7\), as

\[
H=Z^2+\kappa B_3Z+C_6,
\quad B_3=X_0^3+X_0Y_1Y_2+X_1Y_0Y_1.
\]

Here \(C_6\) is determined by the unchanged target and contains the five
ordinary sixth powers and the other three mixed terms. Thus the quotient
is also a double cover of \(\mathbb P^4\), by completing the square.

Let \(G=\operatorname{Gr}(2,5)\), let \(U\) be its tautological rank-two
subbundle, and let \(L=\det U^*\) be the Pluecker line bundle. A weighted
line is a line \(\ell\) in \(\mathbb P^4\) with a cubic
\(Z\in H^0(\ell,\mathcal O_\ell(3))\) satisfying \(H|_\ell=0\).
Consequently its parameter scheme \(M\) is the zero locus of a section
of \(\pi^*\operatorname{Sym}^6 U^*\) on
\(\operatorname{Tot}(\operatorname{Sym}^3 U^*)\). The expected dimension
is \(6+4-7=3\).

This description is independent of a convention for projectivizing a
bundle. The tangent sequence of a vector-bundle total space gives

\[
K_{\operatorname{Tot}(V)}=\pi^*(K_G\otimes\det V^{-1}),
\qquad \det\operatorname{Sym}^m U^*=L^{m(m+1)/2}.
\]

For a regular zero locus, adjunction therefore gives

\[
\omega_M=L^{-5-6+21}|_M=L^{10}|_M.
\]

The morphism \(M\to G\) is finite: locally trivialize the line, evaluate
its cubic at four fixed distinct points, and use the invertible evaluation
matrix. Each value satisfies a monic quadratic with coefficients on the
Grassmannian chart. The four values, hence the cubic's coefficients, are
integral. The remaining incidence equations give a quotient of this
finite algebra. Thus a smooth projective component which is everywhere a
regular zero locus would have ample canonical bundle. The fact that the
ambient quotient is Fano does not make its line space Fano.

The cubic discriminant transforms by \((\det T)^6\) under a change
\(T\) of the two line coordinates, so \(\Delta\) is a section of
\(L^6\). The root model \(\widehat M\), cut out by
\(\delta^6=\Delta\) in \(\operatorname{Tot}(L|_M)\), has the conditional
dualizing-bundle calculation

\[
\omega_{\widehat M}
 =\widehat\pi^*(\omega_M\otimes L^5)
 =\widehat\pi^*L^{15}.
\]

These are adjunction calculations, not a proof that a smooth resolution
of our actual component is of general type. Global regularity,
normality, and discrepancies have not been checked. In particular, the
root model is singular at the seed we actually use.

### The actual root model is singular at the seed

Set \(X_1=1,Y_1=t\) and write

\[
Y_0=\xi+at,\quad X_0=k_0+bt,\quad Y_2=k_2+ct,
\quad Z=k_3+dt+et^2+ft^3.
\]

The seven coefficients of \(H=0\) have invertible Jacobian in
\((\xi,a,b,c,d,e,f)\) at
\((i,0,0,i,0,0,0)\), with \((k_0,k_2,k_3)=(0,0,0)\).
Its determinant in this order is \(-135000i/2401\). Thus \(M\) is smooth
of dimension three at this seed.

Here \(Z=0\), and the discriminant is homogeneous of degree four in
the cubic coefficients. Its differential vanishes. So the new equation
\(\delta^6-\Delta=0\) has zero differential at \(\delta=0\).
Its unnormalized model has tangent dimension four and local dimension
three. It is singular. This is a statement about this model, not about
every model of its function field.

The independently recomputed first incidence jet in endpoint direction
\((\alpha,\beta,\gamma)\) gives

\[
Z=\lambda(\gamma-\alpha t^3)+O(\lambda^2),\qquad
\Delta=-27\alpha^2\gamma^2\lambda^4+O(\lambda^5).
\]

Thus the discriminant is not identically zero on the seed component.
Any argument based on a smooth root cover must deal with this
degeneration explicitly.

## Two alternative quotient screens

For the global parity involution, put
\((Y_0,Y_1,Y_2)=(y,uy,vy)\) and \(r=y^2\). The generic surface fiber
of the quotient over \((u,v)\) is a sextic in
\(\mathbb P(1,1,1,2)\) with coordinates \((X_0,X_1,X_2,r)\).
The exact polynomial is in [the audit output](strategy-geometry.json).

At \(u=v=0\) it is

\[
X_0^6+X_1^6+X_2^6+r^3+\kappa X_0^3X_2^3=0.
\]

The gradient forces \(r=X_1=0\); its remaining equations have no
nonzero solution since \(4-\kappa^2=171/49\ne0\). The surface avoids
the ambient singular point. Thus this fiber is smooth, and the generic
fiber is smooth by openness. Adjunction gives \(K_S=\mathcal O_S(1)\).
Its holomorphic two-forms have residue numerators \(X_0,X_1,X_2\).

Under \(X_2\mapsto\zeta_3X_2\), their characters are respectively
\(\zeta_3,\zeta_3,\zeta_3^2\). The fixed target has character
\(\zeta_3\), so there are **two** holomorphic two-forms in that sector.
This quotient does not meet a proposed zero-geometric-genus criterion
in the target sector. This does not say the particular target restricts
nontrivially or exclude other uses of the quotient.

An exhaustive check of diagonal symmetries recovers the eighteen-element
group from the protected motivation. The target-residue kernel is cyclic
of order six, generated by weights
\((1,5,5,2,2,0)\) modulo six in repository coordinate order. The only
projective involution is global parity: there is no paired-sign
involution of multiplicities two and four.

The invariant linear ratios \(X_1/X_2\) and \(Y_0/Y_1\) give another
natural surface family. Its coordinates have weights \((1,5,2,0)\), and
its invariant holomorphic two-forms include exactly the degree-two
numerators \(X_2^2,Y_1^2\). A smooth specialization at both ratios zero
is checked by the maximum-coordinate gradient bound \(6\kappa<6\).
This second simple vanishing shortcut fails too. The order-three sector
theorem in the book has specific weighted threefold hypotheses; these
counts are not an application of that theorem to our fourfold.

## What a better route must supply

There is a useful geometric interpretation of the current family. The
lift of a quotient line is a plane cubic \(w^3=Z(s,t)\). On its plane,
the sextic factors as

\[
w^6+\kappa B_3w^3+C_6
 =(w^3-Z)(w^3+\kappa B_3+Z),
\]

using the line incidence equation. The generic lifted curve is elliptic,
with the cyclic action and isotrivial geometry already identified by
`benediamond`. The cover of the parameter space is carrying this
elliptic twisting; it is not an auxiliary rational coordinate chart.

The most relevant next candidate is a correspondence using these paired
plane cubics, or a different family not fixed fiberwise by the cubic
deck action. Before a new search, require all of the following:

1. Specify its action on the target character, including any action on
   the parameter space. A deck-invariant family of rational curves with
   trivial parameter action cannot serve the required character through
   its fiber fundamental class.
2. Exhibit an exact seed and prove the expected rank or the needed excess
   dimension. A dimension count by itself is not an existence proof.
3. Compute a nonzero target-sensitive infinitesimal pairing, or an exact
   transport identity for \(A=X_0^3Y_0^3\). Reaching an unrelated sector
   does not suffice.
4. Identify the theorem which turns that calculation into support or a
   rational primitive. Only then assemble a complete target-bound Lean
   proof and run the prescribed instance verifier.

For the present line route, resume coefficient searches only if a pole
analysis or a geometric construction supplies a justified denominator
and a finite ansatz. More unstructured rank failures are not the next
milestone. We have not yet supplied a correspondence passing these gates.

## Reproduction, references, and limits

Run from the repository root:

```sh
python3 research/cyclic-triple-residual/verify_strategy_geometry.py
```

The six exact algebra checks pass using SymPy 1.14.0. They bind the input
record, verify discriminant covariance using generators of GL(2), check
the bundle determinant exponents, recompute the actual seed and its
discriminant jet, and enumerate the quotient characters. The output
distinguishes exact computations from geometric hypotheses. These checks
are not Lean proofs of adjunction or a certificate. Existing Lean sources
are unchanged; the earlier kernel replay is recorded separately in
[lean-checks.json](lean-checks.json).

Prior work and attribution for the target, seed, incidence construction,
and residual thread are listed in [README.md](README.md#fixed-inputs-and-prior-work).
The no-go observation about the deck-invariant rational-curve
correspondence is credited there to `benediamond`; it is not a new claim
of this note. The current subgroup audit reproduces protected curation
data rather than claiming a new symmetry discovery.

For general background on line schemes of double covers, see Appendix A
of Pedro Montero and Ivan Rosas-Soto,
[Algebraic cycles on some Fano varieties. Part I: Quartic double fivefolds,
arXiv:2509.12186v4](https://arxiv.org/abs/2509.12186v4), which credits
Tihomirov and Skauli. Its main theorem concerns a different dimension
and branch degree and is not a theorem for our target. Our canonical
calculation above is derived directly on a vector-bundle total space;
we do not import that appendix's displayed projective-bundle canonical
formula.

This audit provides reasons to change the research strategy. It does not
establish a global birational classification, rule out all rational
currents, or certify the target. No discussion comment or PR was posted.
