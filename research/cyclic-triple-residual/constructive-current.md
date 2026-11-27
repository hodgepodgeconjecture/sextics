# Constructing parameter currents near the seed

Research checkpoint, 2026-11-27, prepared by OpenAI Codex for
`hodgepodgeconjecture`. This supplies an exact rational subtraction and a
formal solution of the remaining scalar equation. It does not supply a
rational solution of the full equation or a target certificate.

## An exact rational subtraction

Retain the fixed source and the notation of [README.md](README.md):

\[
\mathcal D(b)=\sum_i\left(\partial_i b_i-
 \frac{\Delta_i b_i}{6\Delta}\right),\qquad
B=J_s/3-J\Delta_s/(18\Delta).
\]

Choose an endpoint index \(j\) for which \(\Delta_j\ne0\), and set

\[
b_j^{(0)}=\frac{J\Delta_s}{3\Delta_j},\qquad
b_i^{(0)}=0\quad(i\ne j).
\]

Then, exactly in the incidence field,

\[
B-\mathcal D(b^{(0)})=R_j,
\qquad R_j=J_s/3-\partial_j\left(\frac{J\Delta_s}{3\Delta_j}\right).
\]

The explicit \(1/\Delta\) pole is gone. On a smooth discriminant chart
where \(J,\Delta_s\) are regular and \(\Delta_j\) is a unit, this is a
regular residual source. Globally there may be poles on \(\Delta_j=0\)
or in the incidence coordinates; this formula does not remove those.

Writing a correction as \(b=b^{(0)}+\Delta c\), the remaining equation is

\[
\Delta\sum_i\partial_i c_i+\frac56\sum_i\Delta_i c_i=R_j.
\]

This is an exact reformulation with no explicit division by \(\Delta\).
It gives a specified rational current and a specified residual equation,
not a guessed denominator bound or a proof of global solvability.

## A formal radial solution

Work in the smooth projective line-space chart \(k_1=1\), with endpoint
coordinates \(x=(k_0,k_2,k_3)\) centered at the accepted seed. Put
\(R=\sum x_i\partial_{x_i}\). Seek a current \(b_i=x_i f\). Its equation is

\[
R(f)+\left(3-\frac{R(\Delta)}{6\Delta}\right)f=B.
\]

On the blowup, in a generic direction \(x=\lambda(\alpha,\beta,\gamma)\),
write \(\Delta=\lambda^4 D(\lambda)\), where
\(D(0)=-27\alpha^2\gamma^2\ne0\). The equation becomes

\[
\lambda f'+\left(\frac73-\frac{\lambda D'}{6D}\right)f=B.
\]

If \(a(\lambda)=7/3-\lambda D'/(6D)=\sum a_n\lambda^n\), then

\[
f_n=\frac{B_n-\sum_{j=1}^{n}a_j f_{n-j}}{n+7/3}.
\]

In characteristic zero the denominator never vanishes for integral
\(n\ge0\). Thus this constructs a unique formal power-series solution
in the chosen radial gauge. Its coefficients are rational functions of
the generic direction: the incidence coefficients are obtained by the
implicit-function recursion with invertible seed Jacobian. The argument
works in the completion at the generic exceptional divisor, where
\(\alpha\gamma\ne0\); it is not a claim of regularity over the entire
exceptional divisor or at the original three-dimensional origin.

This is different from the earlier unsuccessful homogeneous Euler
integration: that Euler field rescales the affine cone with weights
\((1,1,1,3)\). Here \(k_1\) is fixed and the three endpoint deviations
are scaled equally about a particular smooth point. Its nonresonance
therefore does not contradict the earlier weight-zero obstruction.

## The first current is explicit

Let \(C_3\) denote the rational expression
`compact_C0_lambda3` in [source-deformation-exact.json](source-deformation-exact.json).
The previously proved source expansion is \(B=\lambda^3 C_3/9+O(\lambda^4)\).
The recurrence gives

\[
f=\lambda^3 C_3/48+O(\lambda^4),\qquad
b_i=x_i f.
\]

The script verifies the full three-variable leading identity

\[
\sum_{q\in\{\alpha,\beta,\gamma\}}
\left[\partial_q(qC_3/48)-
 \frac{qC_3\,\partial_q(-27\alpha^2\gamma^2)}
 {288(-27\alpha^2\gamma^2)}\right]=C_3/9.
\]

This selected radial representative has poles in the direction variables
\(\alpha,\gamma\). It illustrates coefficient behavior absent from the
earlier endpoint-polynomial ansatz, but does not prove that every
rational current needs those poles. No finite global ansatz is justified
by this leading term alone.

## Exact reproduction and the remaining gap

```sh
python3 research/cyclic-triple-residual/radial_current.py
```

The script uses the existing target-bound incidence builder. It constructs
24-term incidence and source-variation series over each of
\(\mathbb F_{101}\) and \(\mathbb F_{109}\), in directions \((1,2,3)\)
and \((2,1,4)\). It verifies every incidence and source-variation equation.
Division by the order-four discriminant leaves 20 reliable coefficients;
the radial current equation holds modulo \(\lambda^{20}\) in all four
runs. Every recurrence denominator used is nonzero at these primes.
The leading coefficients agree with the independent characteristic-zero
formula. [radial-current.json](radial-current.json) records the coefficients.
All calculations use exact arithmetic.

The formal existence argument is characteristic zero; the four finite-field
runs check its implementation, not a generic rational identity. The
explicit pole-removal identities and the leading current are checked
symbolically. The original source reduction, seed, and prior authors are
credited in [README.md](README.md#fixed-inputs-and-prior-work).

The next proof obligation is rationality: either construct rational
\(c_i\) in the residual equation above, or identify a rational incidence-field
expression whose expansion is a suitable formal current. Finite jets
cannot settle that obligation. A proposed expression must be substituted
back into the full incidence-field equation before target-bound Lean
assembly. No convergence, rationality, or protected `HasCertificate`
theorem is claimed here. No discussion post or PR was made.
