# Two finite boundary families

Research checkpoint, 2026-11-27, prepared by OpenAI Codex for
`hodgepodgeconjecture`. These exact local calculations complement
[the infinity audit](infinity-boundary.md). They treat two explicit families
of lines lying in \(Y_1=0\); they do not classify every finite boundary
point or prove rational exactness. No Lean certificate is asserted.

Put \(K=k_0,L=k_2,\kappa=5/7,d=\kappa^2-4\). In the anchored chart,
allow the direction of \(Y_1\) to be \(h\), rather than fixing it to 1.
If \(\Delta_{\rm raw}\) is the cubic discriminant in that direction
chart, then

\[
\Delta_{\rm old}=h^{-6}\Delta_{\rm raw},\qquad
\eta=-\frac{h}{3w_0r_{\rm raw}}\,d\xi,\qquad
\widetilde\alpha=\left.\partial_s\eta\right|_\xi,\quad r_{\rm raw}^6=\Delta_{\rm raw}.
\]

Thus \(h=0\) is a boundary of the old slope chart. A zero of the raw
discriminant there need not be a zero of the old discriminant.

## Constant-cubic family

For any \(\rho^6=-1\), there is an exact line

\[
X_0=K,\quad X_1=1,\quad Y_0=\rho(L+t),\quad Y_1=0,
\quad Y_2=L+t,\quad W=w,
\qquad w^2=dK^6-4.
\]

Its cubic \(Z=(w-\kappa K^3)/2\) is constant, so its binary cubic has
a triple zero at infinity. This degeneration must not be treated as an
ordinary simple zero of the discriminant.

Near it write

\[
X_0=K+bt,\quad Y_0=\rho L+x+(\rho+u)t,\quad Y_1=ht,\quad Y_2=L+t.
\]

The Jacobian of \((E_4,E_5,E_6)\) in \((x,u,h)\) at the boundary is

\[
\begin{pmatrix}
-120\rho^5L&-240\rho^5L^2&-4\kappa K^2\rho\\
-24\rho^5&-120\rho^5L&0\\
0&-24\rho^5&0
\end{pmatrix}.
\]

It is invertible for \(K\ne0\). Hence \(b\) is a smooth local parameter.
The weighted initial equations give
\(h=H_4b^4+O(b^5)\), \(x=X_5b^5+O(b^6)\), and
\(u=U_6b^6+O(b^7)\). To specify the coefficients compactly, let
\(W_j\) be the coefficient of \(b^jt^j\) in the cubic Taylor
approximation of \(\sqrt{d(K+bt)^6-4}\), taking endpoint value \(w\):

\[
W_1=3dK^5/w,\quad W_2=(15dK^4-W_1^2)/(2w),\quad
W_3=(20dK^3-2W_1W_2)/(2w).
\]
\[
R_4=15dK^2-2W_1W_3-W_2^2,\quad
R_5=6dK-2W_2W_3,\quad R_6=d-W_3^2.
\]
Then

\[
H_4=\frac{R_4}{4\kappa K^2\rho},\quad
X_5=\frac{R_5-8\kappa K\rho H_4}{24\rho^5},\quad
U_6=\frac{R_6-4\kappa\rho H_4}{24\rho^5}.
\]

The script verifies the weighted equations directly from the fixed target,
and verifies \(H_4,X_5\ne0\) generically by exact coprimality at \(K=1\).
Thus the endpoint map has ramification 5.

The raw cubic coefficients have leading terms

\[
\tfrac12(w-\kappa K^3),\quad
\tfrac12(W_1-3\kappa K^2)b,\quad
\tfrac12(W_2-3\kappa K)b^2,\quad
\tfrac12(W_3-\kappa)b^3.
\]

Their discriminant has order exactly 6. Nonvanishing of its leading
coefficient is checked at \(K=1\), for both roots of
\(w^2+367/49=0\), by polynomial coprimality. Therefore the old
discriminant has order \(6-6\cdot4=-18\). The Kummer cover is unramified
on the normalized branch; \(r_{\rm raw}\) has order 1.

The source deformation adds \(4sX_0^3Y_0^3\) to \(D\). Its forcing
in \((E_4,E_5,E_6)\) has orders at least \((1,2,3)\). Solving at
fixed \(b\) gives

\[
x_s=O(b^2),\qquad u_s=O(b^3),\qquad h_s=O(b).
\]

One can see these bounds without assuming a smooth root-cover family:
scale the unknown derivatives by \((b^2,b^3,b)\) and the three equations
by \((b,b^2,b^3)\). The limiting matrix remains invertible. The necessary
Jacobian-entry orders and source forcing are checked in the script.

Changing to differentiation at fixed \(\xi\) introduces
\(b_s=-x_s/x_b=O(b^{-2})\). It preserves \(h_s=O(b)\), and gives
\(\dot\Delta_{\rm raw}/\Delta_{\rm raw}=O(b^{-3})\): the
discriminant is weighted homogeneous of weight 6 when its cubic
coefficients have weights \((0,1,2,3)\). Since \(d\xi\) has order 4,
\(\eta\) has order 7 and \(\widetilde\alpha\) has order at least 4.
The apparent degeneration produces no pole of the transformed source.

## Translated-cubic family

For \(\rho^6=d/4\), a second exact line is

\[
X_0=K+t,\quad X_1=1,\quad Y_0=\rho(K+t),\quad Y_1=0,
\quad Y_2=L,\quad W=w,\qquad w^2=-4(1+L^6).
\]

Here \(Z=(w-\kappa(K+t)^3)/2\) has discriminant
\(-27\kappa^2w^2/16\ne0\). Use the chart \(b=1\), and write
\(\xi=\rho K+x,a=\rho+u,Y_2=L+ct\). Define

\[
A_5=2\kappa[\kappa L+(\kappa-2)\rho],\qquad
A_4=4KA_5-4\kappa\rho L^2.
\]

The local Jacobian in \((x,u,h)\) is

\[
\begin{pmatrix}
-120\rho^5K&-240\rho^5K^2&A_4\\
-24\rho^5&-120\rho^5K&A_5\\
0&-24\rho^5&0
\end{pmatrix}.
\]

It is generically invertible. Take \(W_1,W_2,W_3\) now from the cubic
Taylor approximation of \(\sqrt{-4(1+(L+ct)^6)}\), and put
\(R_4=-60L^2-2W_1W_3-W_2^2\). The initial equations give

\[
h=H_4c^4+O(c^5),\quad x=X_4c^4+O(c^5),\quad u=U_5c^5+O(c^6),
\]
\[
H_4=\frac{R_4}{KA_5+4\kappa\rho L^2},\qquad
X_4=\frac{A_5H_4}{24\rho^5},\qquad
U_5=\frac{\kappa^2H_4}{12\rho^5}.
\]

At \(K=1,L=2\), exact coprimality with \(\rho^6-d/4\) verifies
the Jacobian and both \(H_4,X_4\) are nonzero at all six values of
\(\rho\). Thus \(c\) is a local parameter and the endpoint
ramification is 4. The old discriminant has order -24, and the
Kummer cover is unramified.

The source forcing is regular. At fixed \(c\), the unit local Jacobian
makes \(x_s,u_s,h_s\) regular. At fixed \(\xi\),
\(c_s=O(c^{-3})\), while \(h_c=O(c^3)\); hence \(h_s\) remains
regular. The raw discriminant logarithmic derivative has pole order at
most 3. Since \(d\xi\) has order 3, \(\eta\) has order 7 and
\(\widetilde\alpha\) has order at least 3. Again it is regular.

## Consequence and remaining gap

For generic \(K,L\), these families supply six ramified boundary
branches of each type, with either endpoint sign. All exclude poles of
a prospective rational primitive at the corresponding points. The
calculations do not assert that these are all lines in \(Y_1=0\).

One possible compact-cohomology test uses the opposite-character
differential \(d\xi/(w_0r_{\rm old}^5)\). It is regular at both families
above and at the infinities already analyzed. At ordinary simple
discriminant zeros it is also regular. Establishing its regularity at
every remaining point, then evaluating the residue-duality obstruction,
would give a substantive test of the one-component primitive. That
global calculation has not been completed. Even an obstruction to this
one-component route would not exclude a three-component current or a
target certificate.

## Reproduction and credit

```sh
python3 research/cyclic-triple-residual/finite_boundary.py
```

All target-binding, boundary, Jacobian, weighted-initial, and
characteristic-zero coprimality checks pass. The formulas and results are
recorded in [finite-boundary.json](finite-boundary.json). The local
normalization and valuation arguments are explained above. Prior
incidence and source constructions are credited in
[README.md](README.md#fixed-inputs-and-prior-work).
No protected file, discussion, or PR is changed by this checkpoint.
