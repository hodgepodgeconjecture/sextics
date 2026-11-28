# A low-degree formulation of the critical trace

Research checkpoint, 2026-11-28, by OpenAI Codex for
`hodgepodgeconjecture`. This is an implementation alternative for the
[principal-parts trace](principal-parts.md), not an evaluation of it.

Fix K,L and k1=1. Let q=(xi,a,b,c,d,e,f,k3), and let Q be the seven
coefficients of 7H. Write J=Q_q and e=(1,0,...,0). On the ordinary chart
where the matrix M with rows J,e is invertible, define the normalized
tangent v by Jv=0, ev=1. The discriminant gradient has a unique decomposition

\[
\nabla\Delta=J^t\lambda+\Delta_\xi e^t.
\]

Thus an ordinary critical point is described by the polynomial equations

\[
Q=0,\qquad J^t\lambda=\nabla\Delta,\qquad Jv=0,\qquad ev=1.
\]

There are 23 variables (eight q, seven lambda, eight v) and 23 equations.
All have total degree at most six after fixing the parameters. Saturation
by det(M), w, Delta, and the second derivative below is still necessary:
the raw equations may include singular or unwanted components. The degree
bound applies to this unsaturated presentation, not to a final Groebner
basis or to additional equations encoding the saturation.

Let A_j be the seven coefficients of the actual source along the line.
Differentiating Q(q(s))-7sA(q(s))=0 at fixed xi gives
J q_s=7A and e q_s=0. Consequently, at a critical point,

\[
\dot\Delta=7\lambda^t A.
\]

In fact the identity holds everywhere if lambda is obtained from the
full gradient decomposition, because e q_s=0. Differentiating Q along xi
twice gives J q''=-Q''[v,v], e q''=0. Therefore

\[
\Delta_{\xi\xi}=v^t\left(\operatorname{Hess}\Delta
 -\sum_j\lambda_j\operatorname{Hess}Q_j\right)v=:B.
\]

The trace summand derived in the preceding note becomes

\[
\frac{14\lambda^t A}{w^2\Delta B}.
\]

This avoids constructing N3, Ns, their expanded determinants, or a
second implicit source derivative. It increases the number of variables;
there is no claim that the resulting elimination is automatically faster.
The tangent and multiplier variables are uniquely determined on the
ordinary chart, so adjoining them does not change the critical algebra
or its multiplication trace there.

## Exact check and limitations

Run `python3 research/cyclic-triple-residual/critical_lagrange.py`.
The program uses the existing incidence builder and the saved modulo-13
source witness. It solves the normalized tangent and multiplier systems,
then independently computes the second derivative through an implicit
acceleration solve. Both methods give Delta_xixi=4, and the multiplier
source formula gives Delta_dot=4. Delta_xi=4 at this witness: it is **not**
a critical point and is not presented as one. This checks the general
identities used to build the critical system. The measured incidence
degrees are all six. Results are in [critical-lagrange.json](critical-lagrange.json).

No complete critical algebra, saturation, or global trace has been
computed. The omitted boundary corrections and global holomorphy of the
test differential remain separate obligations. Even a vanishing corrected
pairing would not prove exactness or settle the three-component current.
Prior incidence and source constructions are credited in
[README.md](README.md#fixed-inputs-and-prior-work). No protected files,
discussion posts, PRs, or cloud compute resources were changed or used.
