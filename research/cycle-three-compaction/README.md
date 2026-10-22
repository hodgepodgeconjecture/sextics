# `cycle-three` (PR #165): what it does, why it is 2.1 MB, and a 4× smaller, faster version

Research notes for `hodgepodgeconjecture`. Nothing here is a solve PR; the
accepted proof in `data/cycle-three/` is unchanged.

## Summary

- **The mathematics of PR #165 is natural and probably close to minimal** for this
  instance. It takes the fibration of the sextic by plane cubic curves, passes to
  their Jacobians, and integrates along rational curves on the resulting rational
  elliptic surfaces. The primitive is then transported back.
- **91% of the 2.14 MB is one classical identity.** That identity is the syzygy
  \(J^2 = 4\Theta^3 + 4s_4\Theta H^4 - s_6H^6 \pmod U\) for the six-coefficient cyclic
  cubic, plus the scaffolding used to check it. The 1.59 MB `CubicIdentityData.lean`
  contains no numerical witness. It writes out, coefficient by coefficient, the
  products that 190 separate `decide` theorems refer to.
- **The same identity can be proved from about 26 KB of Lean in under a minute.**
  The idea is to encode the whole identity as one big-integer equation. Lean's kernel
  does big-number arithmetic with GMP, so it only has to perform a few dozen large
  multiplications. The resulting `full_syzygy` has exactly the original type, so it is
  a drop-in replacement. The same trick replaces the 324 s, 8.75 GB `ring` proof of the
  cubic's differential identity.
- **End to end:** I built a variant `cycle-three` proof with both replacements. It has
  45 modules and 507 KB, against 60 modules and 2.14 MB. It passes the same compile,
  `leanchecker` replay, and audit steps as the protected verifier, with the same final
  theorem and the same three axioms. Measured with identical settings, it compiles in
  917 s instead of 2,188 s and replays in 443 s instead of 1,506 s. Its compile peak
  is 5.1 GB instead of 8.9 GB. Details are in [Results](#results).

---

## 1. The project, for someone new to the Hodge conjecture

**Shapes cut out by equations.** Take six complex variables
\(X_0,X_1,X_2,Y_0,Y_1,Y_2\) and one polynomial \(F\) of degree 6. Its zeros in
projective 5-space form a 4-dimensional complex "shape" \(X\) (8 real dimensions), a
*sextic fourfold*. It is a Calabi–Yau variety: it carries a nowhere-vanishing
holomorphic 4-form \(\omega\).

**Holes and their types.** A donut has 1-dimensional loops. \(X\) has 4-dimensional
"cycles", measured by the cohomology group \(H^4(X)\). Hodge theory sorts each class by
how many holomorphic directions \(dz\) and antiholomorphic directions \(d\bar z\) it uses:
\(H^4 = H^{4,0}\oplus H^{3,1}\oplus H^{2,2}\oplus H^{1,3}\oplus H^{0,4}\). The extreme
piece \(H^{4,0}=\mathbb C\,\omega\) is the most "transcendental" part, the furthest from
geometry.

**The (generalized) Hodge conjecture.** The Hodge conjecture says that rational classes
of type (2,2) come from algebraic subvarieties. Grothendieck's generalized version
covers more. Suppose a rational family of classes avoids the extreme types, so it lives
in \(H^{3,1}\oplus H^{2,2}\oplus H^{1,3}\) ("coniveau 1"). Then it should be
**supported on a divisor**: some hypersurface section \(Y\subset X\) exists such that
every one of these classes becomes zero once \(Y\) is deleted. Informally, the classes
are explained by lower-dimensional geometry.

**This case.** Here \(F\) is even in \(Y\), meaning \(F(X,-Y)=F(X,Y)\). So \(X\) has a
symmetry \(\iota\) that flips the signs of the \(Y\)'s. The classes fixed by \(\iota\)
form a piece \(H\). Since \(\iota\) flips the sign of \(\omega\), \(H\) has no (4,0) part,
so it has coniveau 1 (Voisin). The generalized conjecture therefore predicts a divisor
that supports \(H\). This is open for general \(F\). Voisin showed that it would imply
her conjecture on K3 surfaces that are double covers of the plane.

**From topology to algebra.** By Griffiths and Diamond, the classes in \(H\) are residues
of rational forms \(A\,\Omega/F^{q+1}\), where \(A\) is an odd monomial of degree
\(6q\) and \(q\in\{1,2\}\). "Supported on a divisor \(D\)" is equivalent to finding a
vector field \(v=N/D\) with the following properties:

- \(D\) is not divisible by \(F\);
- each numerator has degree \(\deg D+6q-5\);
- \(F^{q+1}\) divides \(A-q\,dF(v)+F\operatorname{div}(v)\).

Think of it as finding an antiderivative. You subtract an exact form to cancel the pole
along \(X\), and the price is new poles along \(\{D=0\}\).

**The project** fixes many instances \((F,q,A)\). A solve is a Lean 4 proof of that
existence statement. A protected pipeline checks it: a sandboxed compile, a full kernel
replay with `leanchecker`, and an audit of the statement and its axioms. Solving every
odd \(A\) for one \(F\) proves the conjecture for that \(F\).

The difficulty is that nothing bounds \(\deg D\) in advance. Once \(D\) is fixed the
equations are linear, but searching over \(D\) is nonlinear and huge. For `cycle-three`,
`rbajaj5` proved that no certificate has \(\deg D\le 6\).

## 2. What PR #165 does, intuitively

The instance is
\(F=\sum_i (X_i^6+Y_i^6)+\tfrac12\sum_{i \bmod 3} X_i^3X_{i+1}Y_iY_{i+1}\), with
\(q=1\) and \(A=X_0^3Y_0^3\).

1. **Turn the sextic into plane cubic curves.** Pair up the coordinates: put
   \(r_i=Y_i/X_i\) (a "slope") and \(u_i=X_i^2\). Then \(F=0\) becomes a cubic in
   \((u_0:u_1:u_2)\) whose coefficients depend on the slopes:
   \[U=(1+r_0^6)u_0^3+(1+r_1^6)u_1^3+(1+r_2^6)u_2^3+\tfrac12(r_0r_1u_0^2u_1+r_1r_2u_1^2u_2+r_2r_0u_2^2u_0).\]
   Up to a 4-to-1 sign cover, the fourfold is a 3-parameter family of plane cubic curves,
   which are genus-one curves. The target class is \((Y_0/X_0)^3\) times the rate at
   which the curves' periods change when the \(X_0^6\) coefficient moves.
2. **Pass to the Jacobian.** A plane cubic without a chosen point is not an elliptic
   curve in Weierstrass form. Classical invariant theory (Aronhold and Salmon; Fisher's
   normalization) supplies a canonical degree-9 map to its Jacobian
   \(\eta^2=\xi^3+s_4\xi-s_6/4\), with \(\xi=\Theta/H^2\) and \(\eta=-J/(2H^3)\). The
   ingredients are the covariants \(H\) (the Hessian), \(\Theta\), and \(J\). The
   identity behind it is the **syzygy** \(J^2\equiv 4\Theta^3+4s_4\Theta H^4-s_6H^6\)
   modulo \(U\).
3. **The Jacobian side is rational.** A hidden \(\mu_3\)-symmetry makes \(s_4,s_6\) depend
   on the slopes only through \(v^6\). After setting \(T=Rv^6\), the Weierstrass
   coefficients have degrees 4 and 6 in \((R,T)\). For fixed \((c_1,c_2)\), that makes a
   **rational elliptic surface**, and such surfaces are swept out by rational curves.
   The proof integrates along them in three steps:
   - pull a polynomial source back to a family of polynomial curves
     \(T\mapsto(R(T),X(T),Y(T))\);
   - take the antiderivative in \(T\) (a polynomial's antiderivative is a polynomial);
   - push back down by averaging over the finitely many curves through each point, a
     trace.

   The curves come from a 13×13 polynomial system. Instead of solving it, the proof finds
   a nonsingular solution mod 5 and lifts it to characteristic 0 with an étale-lifting
   theorem.

   The original surfaces are not like this. By my calculation, their \(\mu_3\)-quotient is
   a Dolgachev-type surface: \(p_g=q=0\), but Kodaira dimension 1 (\(K\equiv F/3\), two
   triple fibres).
   Such a surface has no families of rational curves, which is why the proof must go
   through the Jacobian. The Jacobian of a Dolgachev surface is a rational elliptic
   surface.
4. **Transport back.** The primitive is pulled back to the sextic. The covariant map
   moves when the coefficient moves, so the derivative picks up a correction: a
   first-variation formula, proved with dual numbers. Two more terms complete the
   transport: a Piola/Jacobian-determinant factor (\(-144\,r\,v^3\)) and a correction for
   changing the normal chart.
5. **Clear denominators** with the shared library. The argument is nonconstructive: no
   explicit \(D\) is produced.

## 3. Where the 2.14 MB goes

| Layer | Files | Bytes | Share |
| --- | ---: | ---: | ---: |
| Cubic covariant syzygy and its encoding machinery | 20 | 1,940,904 | 91% |
| Jacobian family (curves, lifting, primitive) | 16 | 79,263 | 4% |
| Physical transport and clearing | 11 | 71,886 | 3% |
| Cubic map, differential, jets | 9 | 35,626 | 2% |
| Later extracted to the shared library | 3 | 9,692 | 0.5% |

The syzygy is one polynomial identity in \(a,\dots,f,x,y\): \(J^2\) has 22,905 monomials,
the main right-hand side has 21,411, and the \(U\)-multiple has 14,820. Lean's `ring`
cannot handle about 60K terms. The PR's way around this is sound and clever: it uses
Kronecker substitution with weight-graded injectivity and a digit bound at \(2^{42}\).

What made it large is how the products were formed. Every intermediate product was
written out by hand, coefficient by coefficient, as symbolic `Expr` trees: \(H^2,\dots,H^6\),
\(U^2,\dots,U^5\), \(\Theta^2\), \(\Theta^3\), \(J^2\), and 13 quotient terms. That way 190
separate `decide` theorems could each name one coefficient. The same 1.57 MB pattern was
copied into `mixed-biweight` (PR #170) and into the shared library
`HodgeComputer/Constructions/TransitiveCubic`.

## 4. Experiments: three ways to check the same identity

All runs used Lean 4.33.0 with the pinned Mathlib, on 4 cores with 15 GB RAM. Peak
memory is per process and includes about 3.1 GB for importing Mathlib.

| Design | Source for the syzygy layer | Compile | Peak memory | Outcome |
| --- | ---: | ---: | ---: | --- |
| PR #165: materialized products, 190 `decide` theorems | 1.65 MB (+ 30 KB shared) | 506 s wall on 4 cores (1,257 CPU-s) | 5.5 GB | passes |
| A: products formed by the kernel, one declaration | 15 KB | killed after 635 s | > 12 GB | fails (kernel caches every symbolic tree) |
| B: integer summaries (value, masses, grade) per coefficient, one declaration | 18 KB | 777 s (713 CPU-s) | 12.0 GB | passes; too close to the 14 GB budget |
| **C: whole identity as one big-integer equation** | **26 KB** (library 15 KB + application 11 KB) | **32–42 s** | **5.1 GB** | **passes; drop-in `full_syzygy`** |

How design C ([`lean/NineKronecker.lean`](lean/NineKronecker.lean),
[`lean/FullSyzygy.lean`](lean/FullSyzygy.lean)) works:

- **Encoding.** Keep the homogeneous coordinate \(z\). Under the weights
  \(a:(3,0,0,0)\), \(d:(2,1,0,0)\), …, \(x:(0,1,1,1)\), \(y:(1,0,1,1)\), \(z:(1,1,0,1)\),
  everything is homogeneous, and the whole identity has the single grade (36,36,36,18).
  At that grade the exponents of \(a,b,c,z\) are determined by those of \(d,e,f,x,y\),
  which are at most 18. Substituting
  \(a{=}b{=}c{=}z{=}1\), \(d,e,f,x,y = B, B^{19}, B^{361}, B^{6859}, B^{130321}\) is
  therefore injective once \(B=2^{50}\) exceeds the coefficient mass (< \(2^{48}\)).
- **Proofs.** The soundness lemmas (grade checker, injectivity, digit bound) adapt the
  PR's own `CubicKronecker`, `CubicGradeCheck`, and `CubicEncodingCheck` from 6 to 9
  variables. They reuse `PolynomialDigitCheck` and `SignedExpressions` unchanged.
- **The kernel's work.** One `decide +kernel` evaluates both sides at the point:
  roughly 60 GMP multiplications of integers up to about 15 MB each. \(U,H,\Theta,J\) are
  written in homogeneous Horner form so that no huge power is materialized more than
  once.
- **Two pitfalls hit along the way.**
  - A single `omega` call over all the weight equations ran out of memory. The fix was
    small helper lemmas.
  - Concrete exponents such as `X^130321` inside proofs made the kernel unfold
    `Monoid.npow` until it hit "deep recursion". The fix keeps the strides symbolic
    (`X^(strides i)`) in proofs; only the final computation sees numbers.
- **Generator.** [`scripts/generate_fullsyzygy.py`](scripts/generate_fullsyzygy.py)
  regenerates `FullSyzygy.lean` byte for byte from `data/cycle-three/CubicCovariants.lean`.

**The second bottleneck, which the same trick removes.** With the syzygy fixed, the
heaviest module became `CubicDifferential.lean`. It holds one `ring` call proving the
relative identity \(H\{U,\Theta\}-2\Theta\{U,H\}+3J=U\{H,\Theta\}\) (with
\(\{A,B\}=A_xB_y-A_yB_x\)), and it takes 324 s and 8.75 GB. The replacement
[`lean/CubicDifferential.lean`](lean/CubicDifferential.lean) keeps every definition and
the exact statement of `relative_identity`. It checks the homogeneous form
\(\dots+3zJ=\dots\) with one `certify` call. It then instantiates the certificate in
the ring `MvPolynomial (Fin 2) K`, mapping \(x,y,z\) to \(X_0,X_1,1\) and each
coefficient variable \(a\) to \(C\,a\). The result is 110 s and 5.1 GB, and nearly all
of that time is the small `ring` bridge lemmas, not the certificate.

## Results

End-to-end check of a variant `cycle-three` folder. It is built from the accepted files
in three steps. First, the new `FullSyzygy` replaces the old syzygy layer; its only
consumer, `PhysicalCubicJets.lean`, changes one import and one theorem name. Second, in
v2 only, the new `CubicDifferential.lean` replaces the old one. Third, modules no longer
reachable from `Solution.lean` are pruned. The run used
[`scripts/verify_local.py`](scripts/verify_local.py), which repeats
`scripts/proposals/sandbox.verify` step for step: the same `lean -j 2 -M 14336` compile
of every module in import order, a fresh `leanchecker --verbose` per module, and the
trusted `Audit.lean`. Docker isolation was left out because no Docker daemon was
available.

| | Original PR #165 | Variant v1 (new syzygy) | Variant v2 (+ new differential identity) |
| --- | ---: | ---: | ---: |
| Lean modules / source | 60 / 2,137,959 B | 45 / 495,403 B | **45 / 507,151 B** |
| Compile, sum over modules | 2,188 s | 1,164 s\* | **917 s** |
| Compile peak (largest process) | 8.92 GB (`CubicDifferential`) | 8.75 GB (`CubicDifferential`) | **5.12 GB** (`FullSyzygy`) |
| `leanchecker` replay, sum over modules | 1,506 s | 500 s\* | **443 s** |
| Replay peak | 5.20 GB | 5.13 GB | 4.70 GB |
| Trusted audit | `HodgeSubmission.checked`; `propext`, `Quot.sound`, `Classical.choice` | same | same |

\*v1 overlapped with other compiles on the machine. The original and v2 were measured
alone, one after the other, with identical settings. The audit step loads the whole
environment and peaks near 9 GB in every run (original 70 s, v2 52 s); that cost does
not depend on the proof.

Overall, v2 has **4.2× less source**, **compiles 2.4× faster**, **replays 3.4× faster**,
and has a **43% lower compile peak**.

The heaviest modules shift as follows:

- **Original.** `CubicDifferential` compiles in 317 s. The coefficient files C, D, and E
  take about 175 s each to compile, and about 180 s each again in replay.
  `CubicIdentityData` takes 136 s.
- **v2.** No module takes more than 110 s: `CubicDifferential` (110 s, mostly `ring`
  bridge lemmas), `CubicFieldDifferential` (93 s), `PhysicalJacobianCoordinates`
  (67 s), `JacobianSeedBinding` (65 s). The largest replay is 33 s (`FullSyzygy`).

Raw logs and summaries are in [`results/`](results/).

**Reproduce:** `scripts/make_variant.sh <repo> <dir>`, then
`python3 scripts/verify_local.py cycle-three <dir> <copy of formal/lean with the Mathlib cache> <out>`.

## 5. Other ideas I evaluated

| Idea | Verdict |
| --- | --- |
| Verify the syzygy for the physical coefficients \(a=1+r_0^6,\dots\) only | 36× **larger** (\(J^2\): 822,237 terms against 22,905). Proving the identity for generic coefficients, as the PR does, is right. |
| One general ternary-cubic library (10 coefficients), reusable by every cubic family | \(J\) has 83,744 terms against 1,095. Neither PR #165's encoding nor design C scales to it; it would need an equivariance (normal-form) argument instead. |
| Normalize \(d=e=f=1\) with the torus action | Does not help. The weight grading already keeps monomials from merging, so the term count is unchanged. |
| Define \(J\) as the Jacobian determinant of \((U,H,\Theta)\) and use the "Weierstrass ODE" \((\partial\xi)^2=36(\xi^3+s_4\xi-s_6/4)\) on \(U=0\) | Lowers degrees, but needs an extra constant-of-derivation argument. Not worth it now that the identity costs about 40 s. |
| Weierstrass model via a point over a cubic extension (degree-1 map instead of the degree-9 covariant map) | Smaller identities, but the transport layer would have to differentiate the adjoined root. No net gain. |
| Explicit certificate by linear algebra on a guessed support (at \(q=1\), the reduced support \(D_{\rm red}\) suffices) | The construction's support is not explicit: Hessian, trace, and norm divisors, of degree ≳ 24 in physical coordinates. Systems would have 10⁵–10⁶ unknowns, and \(\deg D\ge7\) is already known. A research direction, not a shortcut. |
| Macaulay2 / Singular / Sage | Useful for discovery, but nothing computed outside Lean counts toward acceptance. What helps is choosing certificates the kernel can check cheaply, such as big-integer arithmetic. See [section 8](#8-can-a-cas-discover-a-better-certificate) for experiments. |

## 6. Bottlenecks

1. **No fast, kernel-checkable polynomial arithmetic.** `ring` and `grind` stop scaling
   at about 10⁴ terms, and `native_decide` is banned. The kernel's GMP-accelerated `Nat`
   operations are the fast path. Symbolic trees and sparse lists are slow, and the kernel
   keeps every intermediate value for the whole declaration, so memory runs out (designs
   A and B).
2. **Kernel pitfalls:** concrete huge exponents inside proofs, and large `omega`
   contexts.
3. **CI limits:** 14 GB of memory and 60 minutes. `leanchecker` replay repeats every
   `decide`.
4. **Intrinsic size:** the degree-9 covariant map. The covariant formulas blow up for
   general cubics.
5. **No explicit witness:** nonconstructive proofs produce no \(D\), and extracting one is
   a separate, hard problem.

## 7. Recommendations

- **Shared library (maintenance PR):** add the big-integer identity certificate as a
  module, generic in the weights. Then re-prove `TransitiveCubic`'s syzygy with it, which
  removes 1.57 MB from the shared library.
- **Use the technique in new cubic-fibration solves.** A revision of the accepted
  `cycle-three` proof needs separate owner review, so ask before proposing one.
- **Further trimming of the variant** (507 KB):
  - `CubicReflectionData.lean` (193 KB) and `CubicCovariants.lean` (71 KB) hold two
    copies of the same formulas; one compact copy would do.
  - Import the shared `ResidueTransport` modules instead of the local copies.
- `cycle-q2` (pole order 2) has a different residue character and needs a new
  mechanism. This work does not cover it.

## 8. Can a CAS discover a better certificate?

I tested the three places where Singular or msolve could plausibly replace a
nonconstructive or expensive step. Scripts are in [`cas/`](cas/); they need python-flint,
msolve, and Singular.

**A. An explicit curve family on the Jacobian side, to drop the mod-5 seed, lifting, and
trace.**
- For fixed \((c_1,c_2)\), each Jacobian surface is a degree-6 hypersurface in
  \(\mathbb P(1,1,2,3)\), a del Pezzo surface of degree 1 (its anticanonical pencil is
  the elliptic fibration).
- The degree-\((2,4,6)\) curves are its anticanonical conics (\(-K\cdot C=2\)), and a
  degree-1 del Pezzo surface has 2,160 conic classes. So thousands of curves pass through
  a general base point, and the \(T^0\) equation already forces a square root
  \(y_0^2=X_0^3-\tfrac{27}{4}R_0^6\).
- A full Gröbner solve of the 13×13 incidence system did not finish in 20 minutes
  (Singular, mod 32003).
- The decisive test uses the 240 minimal sections (\(X\) of degree 2, \(Y\) of degree 3
  in \(u\)). Their eliminant over \(\mathbb Q\) is **irreducible of degree 240** at
  \((c_1,c_2)=(2,3)\) and at \((5,7)\) ([`gen_sections.py`](cas/gen_sections.py),
  [`factor_elim.py`](cas/factor_elim.py)). Specialization can only refine a
  factorization, so the 240 sections form **one Galois orbit over
  \(\mathbb Q(c_1,c_2)\)**.
- Frobenius cycle types modulo 60 primes include \(30^8\) (a Coxeter element),
  \(24^{10}\), \(20^{12}\), \(15^{16}\), and \(7^{34}1^2\)
  ([`frobenius.py`](cas/frobenius.py)). This is consistent with the full Weyl group
  \(W(E_8)\).
- Consequences: no section, and hence no conic class, is defined over the parameter field.
  Any explicit curve family lives in an extension of very large degree. The PR's design
  (one nonsingular seed, étale lifting, normalized trace) is the economical route; a CAS
  cannot make it explicit cheaply.

**B. Smaller covariants.** In the chart \(z=1\), \(U\) is linear in \(c\), so
substituting \(c=-(ax^3+by^3+dx^2y+ey^2+fx)\) gives canonical representatives modulo
\(U\). The sizes barely move: \(H\) 14→19, \(\Theta\) 207→274, \(J\) 1,095→1,047.
The syzygy then holds exactly with 12,319 terms in \(J^2\)
([`reduce_mod_u.py`](cas/reduce_mod_u.py)). This is not worth the downstream changes once
the identity costs about 40 s.

**C. An explicit global certificate \((D,N,Q)\).** At \(q=1\), the prescribed-support
theorem means only the support of \(D\) has to be guessed; \(N\) and \(Q\) then follow
by linear algebra. The class to kill, however, lies in the \(E_8\) part of the fibres'
Néron–Severi lattice: it is orthogonal to the zero section and to the fibre. The Galois
group acts on that part like \(W(E_8)\), so a support defined over
\(\mathbb Q(c_1,c_2)\) must, heuristically, contain whole Galois orbits of horizontal
curves (at least 240 per fibre). That would force a very high-degree \(D\). I have not
turned this into a lower-bound proof: one subtlety is that part of a class can also be
killed by exactness in the base directions. Still, it explains two facts: `rbajaj5`'s
bound \(\deg D\ge 7\), and the failure of the earlier Gauss–Manin ansatz with a
degree-42 denominator. Explicit certificates are therefore unlikely to be cheaper to check
in Lean for this instance.

**Where a CAS does help:**
- Designing and validating kernel-friendly certificates: gradings, digit bounds,
  coefficient masses, and Horner forms (section 4).
- Supplying the exact structure of identities: the 13-term syzygy came from Sage's
  formulas.
- Finding cheap seeds: msolve finds mod-\(p\) solutions in well under a second, and a
  seed could be chosen so its inverse matrix is small.
- Producing cofactors for the remaining slow `ring`/`simp` proofs
  (`CubicFieldDifferential`, `PhysicalJacobianCoordinates`, `JacobianSeedBinding`), which
  could then be checked as big-integer certificates.

For other instances, where a small support exists, CAS-driven linear algebra is exactly
how the short explicit solutions in this repository were found.
