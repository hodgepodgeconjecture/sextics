# A discriminant coordinate for the rational-current problem

Research checkpoint, 2026-11-27, by OpenAI Codex for
`hodgepodgeconjecture`. The exact target is unchanged. This identifies
coordinate ramification and an exact differential-field reformulation;
it does not construct a full rational current or a certificate.

## Incidence determinant

Use the seven coefficient equations \(Q=\operatorname{coeff}_t(7H)\)
and dependent variables \(y=(\xi,a,b,c,d,e,f)\). Write
\(M=Q_y\) and \(T=\det M\). The remaining coordinates are the four
endpoints; after fixing \(k_1=1\), the smooth seed component has dimension
three. Implicit endpoint differentiation uses \(M^{-1}\).

The first row of \(M\) is \((42\xi^5,0,0,0,0,0,0)\). Define the
polynomials along the line

\[
P_0=(7H)_{Y_0},\quad P_1=(7H)_{X_0},\quad
P_2=(7H)_{Y_2},\quad W=(7H)_Z.
\]

Then the complementary six-by-six matrix consists of the coefficients
of \(P_0,P_1,P_2,W,tW,t^2W\), in degrees zero through five.
Where \(W\) has degree three with leading coefficient \(w_3\), divide
each \(P_i\) by \(W\). Column subtraction gives

\[
T=42\xi^5 w_3^3\det
 [\operatorname{coeff}_{0,1,2}(P_0\bmod W),
  \operatorname{coeff}_{0,1,2}(P_1\bmod W),
  \operatorname{coeff}_{0,1,2}(P_2\bmod W)].
\]

The script verifies the universal column operation, not just a sampled
determinant equality. The remainder expression is used only on
\(w_3\ne0\); the original six-by-six matrix is valid at the seed as
well, where the degree of \(W\) drops. The identity does not classify
the zero set of \(T\) or assert that the line space is singular there.
It detects failure of this particular endpoint projection to be etale.

## Polynomial cofactor fields

For an endpoint \(k_j\), define a polynomial derivation \(L_j\) on
the incidence coordinate ring by

\[
L_j(y)=-\operatorname{adj}(M)Q_{k_j},\qquad
L_j(k_\ell)=T\,\mathbf1_{j=\ell}.
\]

The adjugate identity gives \(L_j(Q)=0\). On the function field where
\(T\ne0\), this is \(T\partial_{k_j}\).
Let \(A_t\) denote the seven coefficients of the actual source
\(A(\phi)\). For the auxiliary source deformation \(H-sA\), put

\[
S_y=\operatorname{adj}(M)(7A_t),\quad S_k=0,
\qquad N_j=L_j(\Delta),\quad N_s=S(\Delta).
\]

Thus \(\Delta_j=N_j/T\) and \(\Delta_s=N_s/T\). The vector \(S\)
describes source deformation, not a tangent vector to the fixed
incidence variety: \(S(Q)=7T A_t\).

Crucially, the explicit current from the previous checkpoint is

\[
b_j^{(0)}=\frac{J\Delta_s}{3\Delta_j}
          =\frac{J N_s}{3N_j}.
\]

The determinant \(T\) cancels. Adding arbitrary powers of \(T\) to
the previous ansatz is therefore not justified by this current.
Its displayed denominators are \(\xi^5\) (from \(J\)) and \(N_j\).
This is an upper bound for this particular expression, not a necessary
pole divisor for all solutions.

The polynomial \(N_j\) has a concrete geometric meaning: it is the
determinant of the eight-by-eight Jacobian of \((Q,\Delta)\) with
respect to \((y,k_j)\). The block determinant formula gives exactly
the cofactor expression above. It governs the coordinate projection
that replaces \(k_j\) by \(\Delta\). Its zeros need not be intrinsic
singularities of the parameter space.

## Exact one-component equation

Choose \(j=3\), fix \(k_1=1\), and work on \(N_3\ne0\). The
derivation

\[
\delta=\frac{L_3}{N_3}
\]

satisfies \(\delta\Delta=1\), \(\delta k_0=\delta k_2=0\).
It has no displayed denominator \(T\). The equation for a current with
only its third component nonzero is exactly

\[
\delta b-\frac{b}{6\Delta}=\Theta,\qquad
\Theta=\frac{T J_s}{3N_3}-\frac{J N_s}{18N_3\Delta}.
\]

This follows by dividing the original scalar equation by
\(\partial_{k_3}\Delta=N_3/T\); it retains the full source.
The original endpoint Jacobian and \(N_3\) are nonzero in the function
field, as verified on exact Hensel branches. Consequently
\((k_0,k_2,\Delta)\) are generically independent coordinates, and the
incidence field \(E\) is a finite algebraic extension of
\(\mathbb C(k_0,k_2,\Delta)\).

With \(r^6=\Delta\), the equation becomes

\[
\delta(b/r)=\Theta/r.
\]

Thus the one-component route asks for a primitive of the algebraic
differential \(\Theta\,d\Delta/r\) in the function field
\(E(r)\), over the fixed parameter field \(\mathbb C(k_0,k_2)\),
with the appropriate character to reconstruct \(b\in E\).
This is an exact reformulation, not an assertion that algebraic
differentials have rational primitives. The field presentation,
ramification, and exactness still need analysis. A negative answer for
this route would not exclude currents with all three components.

## Leading residual and checks

The characteristic-zero first source jet provides \(\Delta_4\) and
\((\Delta_s)_7\). At the seed, \(J_0=-i/6\), so the leading term
of the pole current is

\[
(b_3^{(0)})_4=\frac{(-i/6)(\Delta_s)_7}
                         {3\partial_\gamma\Delta_4}.
\]

Its regular residual is

\[
R_3=-\frac{\alpha^3}{108}
 +\frac{8\alpha\beta^2+32i\beta^3}{\kappa^4}
 +\frac{6\beta^3\gamma}{\alpha\kappa^3}.
\]

It is nonzero and has degree one in \(\gamma\). The leading correction

\[
c_4=\tfrac32R_3|_{\gamma=0}\,\gamma
     +\tfrac35(\partial_\gamma R_3)\gamma^2
\]

satisfies \(\partial_\gamma c_4-c_4/(3\gamma)=R_3\) exactly.
This verifies the first normal integration step. It does not repair the
failed finite candidate classes or prove that the higher corrections
form a rational function. In particular, the radial representative's
\(1/\gamma\) pole is not forced on this leading current.

Reproduce with:

```sh
python3 research/cyclic-triple-residual/discriminant_chart.py
```

The four symbolic checks pass. At both primes 101 and 109, independently
reconstructed incidence points verify the tangent and source cofactor
identities, the new coordinate Jacobian, and current cancellation.
Their valuations of \((T,N_3,N_s,\Delta)\) are \((0,3,7,4)\).
[discriminant-chart.json](discriminant-chart.json) records the results.
The prior incidence builder and its fixed-target binding are unchanged;
prior work and attribution are listed in
[README.md](README.md#fixed-inputs-and-prior-work).

No unrestricted denominator bound, global singular-locus classification,
rational primitive, or protected Lean certificate has been obtained.
No discussion post or PR was made.
