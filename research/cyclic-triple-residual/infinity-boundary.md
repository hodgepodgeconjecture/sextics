# Infinity of the one-component integration curve

Research checkpoint, 2026-11-27, prepared by OpenAI Codex for
`hodgepodgeconjecture`. This continues [the explicit curve model](algebraic-curve-model.md).
It proves local regularity of the transformed source above endpoint infinity
for generic endpoint parameters. It does not prove rational exactness, classify
the finite boundary, or supply a target certificate. The argument is algebraic;
it has not been formalized in Lean.

## Chart and the excess family

Put \(K=k_0,L=k_2,\kappa=5/7,\epsilon=1/\xi\). Work in the chart \(Y_0=1\):

\[
X_0=K\epsilon+bt,\quad X_1=\epsilon+vt,\quad
Y_1=ht,\quad Y_2=L\epsilon+ct.
\]

The direction is \([b:v:c:h]\). Write \(D=\kappa^2B^2-4C\), with
the same fixed \(B,C\) as in the curve model, and reconstruct \(W\) from
\(W^2=D\). Its endpoint value is a unit, tending to \(2i\) at one infinity
and \(-2i\) at the other. Set
\(E_4=D_4-2W_1W_3-W_2^2\),
\(E_5=D_5-2W_2W_3\), and \(E_6=D_6-W_3^2\).
At \(\epsilon=0\), all three lower \(D_j\), \(1\le j\le3\), vanish,
so these equations reduce to \(D_4,D_5,D_6\).

There is an excess plane sextic at \(h=0\):

\[
p_6=(\kappa^2-4)b^6-4(v^6+c^6)=0.
\]

Resolve horizontal branches by writing \(h=\epsilon H\). The initial
equations are \(E_4/\epsilon^2=p_4+O(\epsilon)\),
\(E_5/\epsilon=p_5+O(\epsilon)\), \(E_6=p_6+O(\epsilon)\), where

\[
p_4=A H^2+QH+R,\qquad p_5=UH+V,
\]
\[
\begin{aligned}
A&=\kappa^2v^2,\\
Q&=6\kappa(\kappa-2)Kb^2v+2\kappa(\kappa-2)b^3
       -8\kappa Kbc^2-8\kappa Lb^2c,\\
R&=15(\kappa^2-4)K^2b^4-60v^4-60L^2c^4,\\
U&=2\kappa(\kappa-2)b^3v-4\kappa b^2c^2,\\
V&=6(\kappa^2-4)Kb^5-24v^5-24Lc^5.
\end{aligned}
\]

Eliminating \(H\) gives the degree-twelve selector
\(S=AV^2-QUV+RU^2\). At \(K=1,L=2\), exact computations modulo
1009 and 1013 prove the following for \(p_6=S=0\):

- The affine quotient in \(c=1\) has dimension 72. Its degree-16 and
  degree-17 relation matrices have ranks 81 and 99, respectively.
- The multiplication matrices commute, satisfy the equations, and represent
  the advertised monomial basis.
- The selector Jacobian and \(b,v,U,V\) are units in that quotient.
- There are no points at \(c=0\); the degree-17 binary relation matrix
  has full rank 18.

The nonzero minors give the same conditions on a nonempty characteristic-zero
open set of \((K,L)\). Bezout and the Jacobian criterion give 72 distinct
points, each with \(H=-V/U\ne0\). The implicit function theorem gives 72
unramified branches in \(\epsilon\).

## The other 48 sheets

In the chart \(h=1\), the limiting equations are

\[
\begin{aligned}
D_4&=\kappa^2v^2,\\
D_5&=2\kappa(\kappa-2)b^3v+2\kappa^2bcv-4\kappa b^2c^2,\\
D_6&=\kappa^2b^2(b^2+c)^2
 -4(b^6+v^6+1+c^6+\kappa bvc).
\end{aligned}
\]

Thus there are twelve isolated limiting points. At six of them,
\(v=c=0\) and \((\kappa^2-4)b^6=4\). Eliminate \(b\) by the unit
derivative of \(E_6\). The weighted initial equations, with weights
\((\epsilon,v,c)=(4,2,1)\), are

\[
2\kappa(\kappa-2)b^3v-4\kappa b^2c^2=0,\qquad
\kappa^2v^2+2\kappa(\kappa-2)b^3\epsilon=0.
\]

Here \(b\) denotes its nonzero limiting value. They give one smooth branch
with parameter \(c\), with \(v\) of order 2 and \(\epsilon\) of order 4.
The slope Jacobian at fixed \(\epsilon\) has order 3.

At the other six points, \(v=b=0,c^6=-1\). Eliminate \(c\), again by
a unit derivative of \(E_6\). With weights
\((\epsilon,b,v)=(4,2,3)\), the initial equations are

\[
-4\kappa c^2b^2-24Lc^5\epsilon=0,\qquad
\kappa^2v^2+2\kappa c(\kappa-4Kc)b\epsilon=0.
\]

For generic \(K,L\), eliminating \(\epsilon\) gives a cusp
\(v^2=\text{unit}\cdot b^3+\text{higher weighted terms}\). Its single
normalized branch has orders \((\epsilon,b,v)=(4,2,3)\); its slope
Jacobian has order 5. The required nonzero coefficient condition holds
at \(K=1,L=2\), since \(c^6+1\) and \(5/7-4c\) are coprime.

These initial forms are checked against the full \(D_j\). The corrections
from the reconstructed \(W_j\) have strictly higher weight. The checks
retain the \(2\kappa^2bcv\) term in \(D_5\) and the
\(2\kappa^2bc\epsilon\) term in \(D_4\); dropping these terms would
give an incorrect displayed coefficient in the second balance.

These twelve degree-four branches account for another 48 sheets. Together
with the 72 unramified branches they exhaust the already established
generic degree 120. This degree count also rules out further horizontal
branches over this infinity. The same analysis applies at the other sign
of the endpoint square root.

## Source valuations

Let \(\Delta_{\rm raw}\) be the discriminant of the cubic in the new chart,
and \(r_{\rm raw}^6=\Delta_{\rm raw}\). Changing the binary line coordinates
back to the anchored chart gives

\[
\Delta_{\rm old}=(\epsilon h)^{-6}\Delta_{\rm raw},\qquad
r_{\rm old}=r_{\rm raw}/(\epsilon h).
\]

The transformed source is the source derivative, at fixed endpoint
coordinate, of the differential

\[
\eta=-\frac{d\xi}{3w_{\rm old}r_{\rm old}}
 =\frac{\epsilon^2 h}{3w_{\rm raw}r_{\rm raw}}\,d\epsilon.
\]

The equality \(\widetilde\alpha=\partial_s\eta\) follows directly by
\(w_{{\rm old},s}=2K^3\xi^3/w_{\rm old}\). The auxiliary deformation
remains exactly \(H-sX_0^3Y_0^3\). In the new chart it adds
\(4s(K\epsilon+bt)^3\) to \(D\); it does not change the target.

For the 72 branches, \(h=\epsilon H\) and \(b,H\) are units.
The orders of \(W_1,W_2,W_3\) are at least \(5,4,3\), and their
explicit source derivatives have orders at least \(2,1,0\).
The source forcing of \((E_4,E_5,E_6)\) has orders at least \((5,4,3)\).
After the blowup divisions all are at least 3. The invertible reduced
Jacobian therefore makes all slope derivatives regular, of order at
least 3. At the limit,
\(\Delta_{\rm raw}=27\kappa^2b^6/4\ne0\).
Thus \(\eta\) and \(\widetilde\alpha\) have order at least 3.
In fact the latter's leading coefficient at the \(2i\) endpoint is
\(-H/(18\kappa r_{\rm raw})\ne0\).

At the first isolated type, \(\Delta_{\rm raw}\) is a unit and the
slope Jacobian has order 3. All source forcing and Jacobian cofactors are
regular, so all slope derivatives have poles of order at most 3.
The order of \(\eta\) is 11 and that of its source derivative is at
least 8.

At the cusp type, the cubic coefficients have orders \((0,4,3,2)\).
Consequently \(\Delta_{\rm raw}\) has order 4, with the unique leading
term \(-27Z_3^2Z_0^2\). The slope derivatives have poles of order at
most 5. Every discriminant partial derivative has order at least 2;
hence \(\dot\Delta_{\rm raw}/\Delta_{\rm raw}\) has pole order at
most 7. Before the Kummer extension the order of \(\eta\) is \(31/3\),
and the source order is at least \(10/3\). The Kummer ramification is
3; after writing \(u=z^3\), the source differential has order at least
\(3(10/3)+2=12\).

| Type | Branches | Endpoint ramification | Order of old discriminant | Source order on the Kummer normalization |
| --- | ---: | ---: | ---: | ---: |
| Excess-family selection | 72 | 1 | -12 | at least 3 |
| \(c=0,b\ne0\) | 6 | 4 | -24 | at least 8 |
| \(b=0,c\ne0\) | 6 | 4 | -20 | at least 12 |

All these source differentials are regular. A rational primitive, if it
exists, cannot have a pole at any of these points. This is a useful boundary
condition for a complete rational integration calculation, not an existence
argument. Regularity at infinity and zero residues leave compact de Rham
cohomology unresolved.

## Reproduction and attribution

```sh
python3 research/cyclic-triple-residual/infinity_boundary.py
```

The symbolic leading forms, correction-order inequalities, and two finite
algebra certificates pass. [The JSON report](infinity-boundary.json) records
their results. The normalization and source-valuation arguments are given
above, rather than represented as a machine-checked global theorem.
Prior work on the quotient and incidence construction is credited in
[README.md](README.md#fixed-inputs-and-prior-work). No protected target
files, discussion posts, or PRs are changed by this checkpoint.
