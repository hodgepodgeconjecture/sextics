# Paired cubics: trace and connection obstructions

Research checkpoint, 2026-11-27. Prepared by OpenAI Codex for
`hodgepodgeconjecture`. No target certificate or unrestricted
nonexistence result is claimed. This tests the specific candidate proposed
in [the strategy reassessment](strategy-reassessment.md).

## Two precise conclusions

The natural exchange correspondence through the cubic quotient annihilates
the target character. A second possibility, identifying the compact
elliptic connections of the two cubics over the line-parameter field,
also fails: an exceptional-divisor valuation would force a nonintegral
order for a rational coefficient. Neither argument excludes arbitrary
parameter currents, other correspondences, or a certificate for the fixed
target.

## The natural exchange loses the character

Write the fixed quotient equation as

\[
H(Z)=Z^2+\kappa B_3Z+C_6,\qquad \kappa=5/7.
\]

Its involution sends \(Z\) to \(Z'=-\kappa B_3-Z\). On a quotient
line satisfying \(H(Z)=0\), the plane section factors as

\[
H(w^3)=(w^3-Z)(w^3-Z').
\]

Let \(q:X\to Y\) be the cubic quotient and let \(\tau:Y\to Y\)
be the involution. The natural multivalued correspondence sending a point
to the three points above its exchanged quotient point acts by
\(q^*\tau^*q_*\). The source residue has character \(\zeta_3\), so
its trace under \(q\) vanishes. Indeed, averaging gives the factor
\(1+\zeta_3+\zeta_3^2=0\).

The same cancellation is visible in the transported differential:
both \(w^{-2}\) and \(w^{-5}\) have character \(\zeta_3\).
Thus neither summing the three lifts nor choosing an algebraic lift and
then taking that same full trace preserves the desired source. On each
pair of elliptic curves, the fiber-product correspondence over their
common \(\mathbb P^1\) likewise factors through the trace to the line.
This is not a nonvanishing correspondence in the required sector.

## A local valuation prevents rational identification

Use the actual smooth line-space chart \(X_1=1,Y_1=t\) from the previous
audit. Its seed is
\((Y_0,X_0,Y_2,Z)=(i,0,it,0)\). The seven incidence equations have
invertible Jacobian in the seven dependent coefficients. In endpoint
direction \((\alpha,\beta,\gamma)\), the exact first implicit jet is

\[
\begin{aligned}
Y_0&=i-\kappa\gamma\lambda t/6+O(\lambda^2),\\
X_0&=\lambda(\alpha-6\beta t/\kappa)+O(\lambda^2),\\
Y_2&=it+\beta\lambda+O(\lambda^2),\\
Z_+&=\lambda(\gamma-\alpha t^3)+O(\lambda^2).
\end{aligned}
\]

Substitution into \(Z_-=-\kappa B_3-Z_+\) gives

\[
Z_-=-i\kappa t+\lambda\left[-\gamma+
 (-i\kappa\alpha+\kappa^2\gamma/6)t^2+
 (\alpha+6i\beta)t^3\right]+O(\lambda^2).
\]

For the binary cubic discriminants,

\[
\begin{aligned}
\Delta_+&=-27\alpha^2\gamma^2\lambda^4+O(\lambda^5),\\
\Delta_-&=-4i\kappa^3(\alpha+6i\beta)\lambda+O(\lambda^2).
\end{aligned}
\]

These leading terms are unaffected by the omitted higher incidence jets:
the discriminant is quartic at the zero cubic, while its differential at
\(-i\kappa t\) gives the displayed linear term for the partner.

Take \(\alpha\gamma(\alpha+6i\beta)\ne0\) generically. Blowing up
the smooth three-dimensional line space at the seed gives a divisorial
valuation \(v\) of its function field \(E\), with
\(v(\Delta_+)=4\), \(v(\Delta_-)=1\). This uses a generic direction
on the exceptional divisor, not a numerical arc that could lie in the
zero set of an unknown denominator. Every element of \(E^*\) has
integer valuation.

The compact rank-one connections already derived in this research are

\[
\nabla\omega_\pm=-\tfrac16d\log\Delta_\pm\otimes\omega_\pm.
\]

For a horizontal map \(\omega_+\mapsto h\omega_-\), with
\(h\in E^*\), the coefficient must satisfy

\[
d\log h=\tfrac16(d\log\Delta_- -d\log\Delta_+).
\]

Taking the logarithmic residue along the exceptional divisor requires
\(v(h)=(1-4)/6=-1/2\), a contradiction. Equivalently,
\(\Delta_-/\Delta_+\) cannot be a sixth power times a nonzero constant
in \(E\).

The opposite compact eigenspace is dual to this line under the flat
elliptic polarization. Mapping into that dual would instead require
\(v(h)=-(1+4)/6=-5/6\). Thus mixing the two compact deck eigenspaces
does not repair a rational horizontal map. This is an obstruction to
these maps of connections, not to all algebraic operations on the total
incidence space.

After adjoining a sixth root of \(\Delta_-/\Delta_+\), a horizontal
coefficient can exist. Such an extension must have even ramification
index at this valuation; for the dual identification the index must be
divisible by six. This is compatible with the elliptic curves becoming
isomorphic over a larger algebraic field. It does not provide a rational
current on \(E\), nor a nonzero descended correspondence. In particular,
a traced identification which descended to a horizontal map of these
connections over \(E\) would have to be zero by the argument above.

## Consequence for the research plan

The paired-cubic factorization alone supplies no shortcut to the target.
We should not spend further effort formalizing an exchange argument that
loses the source character, or guessing a rational gauge that cannot
exist. A replacement must either use genuine parameter-direction
information beyond a fiberwise identification, or use a different
correspondence. The scalar-current equation remains open.

The underlying quotient, seed, and isotrivial elliptic interpretation
come from the prior work credited in
[README.md](README.md#fixed-inputs-and-prior-work), including `rbajaj5`,
`fqlx`, and `benediamond`. The universal connection is recorded in
[connection.json](connection.json) and checked by
`verify_cubic_connection.py`. The calculations here are a direct test of
the candidate on that fixed family; no novelty claim is made about
unreviewed external attachments.

Reproduce from the repository root:

```sh
python3 research/cyclic-triple-residual/verify_paired_cubics.py
```

Five exact algebra checks pass; their output is
[paired-cubics.json](paired-cubics.json). The script recomputes the
implicit first jet from the exact target and checks why its truncation
suffices for the discriminant valuations. The valuation and connection
arguments above are mathematical deductions from those identities,
not new Lean-certified theorems. No protected file was edited; no
instance-verifier run, discussion post, or PR was made.
