# An explicit curve model and a simpler source differential

Research checkpoint, 2026-11-27, by OpenAI Codex for
`hodgepodgeconjecture`. This tackles the algebraic integration formulation
by supplying its field equations and transporting the exact source.
No primitive or protected target certificate is claimed.

## Endpoint curve

Fix \(k_1=1\), regard \(k_0,k_2\) as parameters, and put
\(K=\mathbb C(k_0,k_2)\). The constant incidence coefficient is

\[
\xi^6+k_0^6+1+k_2^6+k_3^2+\kappa k_0^3k_3=0,
\qquad\kappa=5/7.
\]

With \(w_0=2k_3+\kappa k_0^3\), this is

\[
w_0^2=-4(\xi^6+C),\qquad
C=1+k_2^6+(1-\kappa^2/4)k_0^6.
\]

Generically \(C\ne0\); the six finite branch points of this double
cover of the \(\xi\)-line are simple. Its smooth compactification has
genus two by Riemann--Hurwitz. This is the endpoint curve, not yet the
full line-incidence curve or its Kummer cover.

## Eliminate three lift coefficients

Write \(H=Z^2+\kappa B_3 Z+C_6\) as before. On a line use

\[
Y_0=\xi+at,\quad X_0=k_0+bt,\quad Y_1=ht,
\quad X_1=1,\quad Y_2=k_2+ct.
\]

The homogeneous direction coordinates are \([a:b:c:h]\); the working
chart is \(h=1\). Set

\[
D(t)=\kappa^2B_3(t)^2-4C_6(t)=\sum_{i=0}^6D_i t^i,
\qquad W(t)=2Z(t)+\kappa B_3(t).
\]

On \(w_0\ne0\), the equation \(W^2=D\), with \(W(0)=w_0\),
determines the three other coefficients of \(W\):

\[
n_2=4w_0^2D_2-D_1^2,\qquad
n_3=8w_0^4D_3-4w_0^2D_1D_2+D_1^3,
\]
\[
W_1=\frac{D_1}{2w_0},\qquad
W_2=\frac{n_2}{8w_0^3},\qquad
W_3=\frac{n_3}{16w_0^5}.
\]

The remaining equations are exactly

\[
\begin{aligned}
P_4&=64w_0^6D_4-4D_1n_3-n_2^2=0,\\
P_5&=64w_0^8D_5-n_2n_3=0,\\
P_6&=256w_0^{10}D_6-n_3^2=0.
\end{aligned}
\]

They are homogeneous of degrees 4, 5, and 6 in \((a,b,c,h)\).
Recover \(Z=(W-\kappa B_3)/2\), hence \(d,e,f\), rationally.
Thus the seven-equation incidence model is replaced, on this open set,
by the endpoint curve and three equations in the three affine slopes.
All \(D_i\), \(B_i\), source coefficients, and residual templates are
exported in [algebraic-curve-model.json](algebraic-curve-model.json).

The accepted seed has \(w_0=0\); this elimination chart does not contain
the seed itself. Its generic nearby incidence points do lie in the chart.
The reconstruction is checked on those points, rather than applying the
fractions at \(w_0=0\).

## Generic degree over the endpoint curve

For each of the two specializations below, the degree-13 multiples of
\(P_4,P_5,P_6\) restricted to \(h=0\) span all homogeneous polynomials
of degree 13 in \(a,b,c\).

| Prime | \((\xi,k_0,k_2,k_3)\) | \(w_0\) | Matrix shape | Rank |
| --- | --- | ---: | --- | ---: |
| 101 | \((1,1,2,13)\) | 70 | \(136\times105\) | 105 |
| 109 | \((1,1,2,58)\) | 70 | \(136\times105\) | 105 |

This proves there are no common projective zeros at \(h=0\) at these
specializations. The report records 105 independent original rows for
replay. These are exact finite-field rank calculations, not approximate
root searches.

Each endpoint satisfies the endpoint equation and has nonzero derivative
\(w_0\) in \(k_3\). It therefore lifts to a characteristic-zero
algebraic endpoint with the same residue. A corresponding nonzero minor
remains nonzero in characteristic zero, and consequently on a generic
open subset of the endpoint parameter space.

Three hypersurfaces of degrees 4, 5, and 6 in \(\mathbb P^3\), with
no intersection at infinity in this affine chart, have a zero-dimensional
complete intersection of length \(4\cdot5\cdot6=120\). Hence the total
generic line-incidence degree over the endpoint curve is 120, counting
all sheets and multiplicities. We have not proved irreducibility: the
function field of the chosen component has degree **at most 120** over
\(K(\xi,w_0)\). Adjoining \(r^6=\Delta\) gives degree at most 720 over
that endpoint field, or at most 1440 over \(K(\xi)\).

These are not degree bounds over \(K(\Delta)\); that projection needs
separate analysis. Nor does a degree bound prove exactness of a
differential.

## Reduced Jacobian and preservation of the source

Let \(S=\partial(P_4,P_5,P_6)/\partial(a,b,c)\), at \(h=1\), and
let \(T\) be the previous seven-by-seven Jacobian for coefficients of
\(7H\). Eliminating \(\xi,d,e,f\) gives, on the incidence equations,

\[
2^{26}w_0^{21}T=42\cdot7^6\xi^5\det S.
\]

The factor and sign follow by taking the Schur complement: the first
four pivots are \(42\xi^5,7w_0,7w_0,7w_0\); permuting the three lift
columns past the slope columns contributes a minus sign, which cancels
the three minus signs from \(7H_i=-7P_i/(4\,\text{denominator}_i)\).
The script checks the cleared identity on independently reconstructed
Hensel branches at both primes to input precision 48. The valuations are
\(v(w_0)=1,v(T)=0,v(\det S)=21\). It also checks reconstruction of
all lift coefficients and all three reduced equations there.

The fixed source is retained: \(H-sA\) replaces \(D(t)\) by
\(D(t)+4sA(t)\), where \(A(t)=X_0(t)^3Y_0(t)^3\). This is verified
as a polynomial identity. At fixed \(\xi,k_0,k_2\), the endpoint
variation is \(\dot w_0=2k_0^3\xi^3/w_0\). The three slope variations
can now be obtained from

\[
S\begin{pmatrix}\dot a\\\dot b\\\dot c\end{pmatrix}
=-\begin{pmatrix}
\partial_{w_0}P_4\,\dot w_0+4\sum_{i=1}^6\partial_{D_i}P_4 A_i\\
\partial_{w_0}P_5\,\dot w_0+4\sum_{i=1}^6\partial_{D_i}P_5 A_i\\
\partial_{w_0}P_6\,\dot w_0+4\sum_{i=1}^6\partial_{D_i}P_6 A_i
\end{pmatrix}.
\]

In this expression the \(D_i\) are treated as independent template
variables for the displayed partial derivatives; slope dependence is in
\(S\). Reconstructing \(Z\) then gives \(\dot\Delta\) rationally.
The dot means source deformation at fixed \(\xi\), unlike the earlier
\(\Delta_s\) at fixed \(k_3\).

## Exact subtraction of the endpoint pole

The differential whose primitive is sought is
\(\alpha=B\,dk_3/r\), with \(r^6=\Delta\). The endpoint relation
gives \(dk_3/d\xi=-6\xi^5/w_0\). Also

\[
\Delta_s=\dot\Delta+\frac{k_0^3}{6\xi^2}\Delta_\xi.
\]

Substitution and exact differentiation give

\[
\alpha+d\left(\frac{k_0^3}{18\xi^2 w_0r}\right)
=\left(\frac{2k_0^3\xi^3}{3w_0^3r}
 +\frac{\dot\Delta}{18w_0\Delta r}\right)d\xi
=:\widetilde\alpha.
\]

This identity is checked symbolically, using
\(w_{0,\xi}=-12\xi^5/w_0\) and
\(r_\xi=r\Delta_\xi/(6\Delta)\). The new density has no displayed
negative power of \(\xi\); the variation formula also introduces none
where \(w_0\), \(\Delta\), and \(\det S\) are units. These conditions
matter: this is not a claim that the differential is globally regular.

If \(G\) is a primitive of \(\widetilde\alpha\) with the required
Kummer character, then the original current is

\[
b=rG-\frac{k_0^3}{18\xi^2w_0}.
\]

The missing step is precisely such a rational primitive in the specified
algebraic function field. The exact subtraction and explicit curve model
do not establish its existence.

## Reproduction and scope

```sh
python3 research/cyclic-triple-residual/algebraic_curve_model.py
```

The symbolic identities, rank certificates, and local component checks
all pass. This gives a concrete model for further ramification and
differential reduction, not a completed integration algorithm or a
Lean-certified theorem. Prior work on the quotient, seed, incidence,
and source transport is credited in
[README.md](README.md#fixed-inputs-and-prior-work). No protected file,
discussion, or PR was changed by this checkpoint.
