# Target research playbook

How we plan to pick and attack the next open target, and how to work differently from
the first solves. It combines three sources:

- the [cycle-three compaction study](cycle-three-compaction/README.md);
- the CAS experiments in [cycle-three-compaction/cas](cycle-three-compaction/cas/);
- a survey of upstream as of 2026-11-26 (book chapters, shared library, open-target notes).

Notes for `hodgepodgeconjecture`. Everything marked **hypothesis** is an idea to test,
not an established result.

---

## 1. Ten principles

1. **Check for existing work first.** Look for open PRs, upstream branches (for
   example `stack/pr221-cycle-q2`), and the instance discussion before computing
   anything.
2. **Build the character ledger before choosing a method.** About half the open
   targets lie on a variety where another sector is already solved. The only
   obstacle there is that the target's character is not the one an existing quotient
   serves.
3. **Classify the fibre geometry before searching for curves.**
   - Rational or uniruled fibres: integrate along rational curves directly.
   - Fibres with \(p_g=0\) but Kodaira dimension 1 (Dolgachev-type): go through the
     Jacobian.
   - Fibres with \(p_g>0\) in the needed sector: that route is closed.
4. **Run the Galois screen before searching for an explicit certificate.** It takes
   seconds. A huge Galois orbit means no small explicit support exists, so use a seed,
   lift it, and take a trace. A small orbit suggests trying explicit linear algebra.
5. **Treat formalizability as a filter, not an afterthought.** The book proves
   existence for several open targets with Hodge-theoretic arguments that Lean cannot
   reach (Abel–Jacobi maps, Hodge loci, the Hodge conjecture on an auxiliary
   variety). For those targets the work is to find an algebraic route.
6. **Prove identities for a whole coefficient family, then specialize.**
   Specializing first made the `cycle-three` syzygy 36× larger.
7. **Memory is the binding Lean limit.** The budget is now 5.5 h, 996 files, and
   14 GB. Design certificates so that each declaration's kernel working set stays
   small.
8. **Prototype the most expensive identity check first (risk first).** If it cannot be
   checked cheaply, change the construction before writing the rest.
9. **Generate large Lean data with scripts that reproduce it byte for byte.** Check the
   generators in next to the proof.
10. **Share partial and negative results in the instance discussion,** with the search
    space stated precisely. The project asks for this, and it stops others repeating
    the work.

---

## 2. The landscape (2026-11-26)

- **Instances:** 82 in total, 67 solved, 15 open.
- **Solve policy** (`.github/solve-policy.json`): at most 996 files, a 19,800 s (5.5 h)
  time budget, and 14,336 MB of memory. The time budget used to be 60 minutes.
- **Toolchain:** Lean 4.36.0 (#267). The files in `research/cycle-three-compaction/lean`
  were checked under 4.33.0 and need re-checking.
- **Shared library highlights** (`formal/lean/HodgeComputer/Constructions/`):
  - `PolynomialEncoding`: signed expressions, bivariate polynomials, the digit check,
    `KroneckerCheck`, and `PolyOps`. The last two give multi-variable mixed-radix
    encodings and the factored Kronecker × grid check of book §5.8.
  - `QuadraticCovers`: pair algebra, Laurent jets, and rational cover residues. This is
    the residue-cancellation machinery for rational sources used in
    `symmetric-sextic-residual`.
  - `PolynomialCurves`: arithmetic lifting and the polynomial-family certificate.
  - `ResidueTransport`: change of density, the moving map, and coordinate independence.
  - `HigherResidues`, `EllipticCompanions`, `HoldoutCover`, `RationalFields`,
    `FiniteAlgebra`, `BinaryCubic`, `CubicSurfaces`.
  - `TransitiveCubic` is still 2.0 MB of per-coefficient tables. The graded
    whole-identity certificate from the `cycle-three` study would shrink it (§6.2).
- **New book material worth reading before any target:**
  - [§3B.11](../docs/03b-weighted-threefolds.md): moving lines and global residue
    cancellation.
  - [Chapter 3D](../docs/03d-hodge-loci.md): Hodge loci, sector counts (§3D.6), and
    paired-exchange families (§3D.8).
  - [§4.6](../docs/04-known-locus-and-frontier.md): extension obstructions and the
    transcendental part of surface fibres.
  - [§5.5](../docs/05-computation.md): rational exactness on curves, as a finite
    Riemann–Roch test.
  - [§5.8](../docs/05-computation.md): Kronecker-grid verification.
- **Existing screens:**
  - [`scripts/audit-jacobian-scope.py`](../scripts/audit-jacobian-scope.py): χ(O) and
    Kodaira dimension of paired-cubic Jacobian models.
  - [`cycle-three-compaction/cas`](cycle-three-compaction/cas/): Galois-orbit and
    Frobenius screens.

---

## 3. Open targets at a glance

Discussion numbers come from `discussions/index.json`. The notes summarize each
target's `why.md` as of 2026-11-26; read the file before starting.

### A. Existence proved in the book, formal proof missing

The mathematics says a certificate exists. What is missing is an algebraic,
formalizable construction.

| Target | Discussion | What is known | Main obstruction | Leads |
| --- | --- | --- | --- | --- |
| `cycle-q2` (q=2) | #61 | §3D.5 (a paired sign plus an order-3 symmetry) supports the target's sector, of type (6,26,6). | Moving curves away from the Fermat seed is non-explicit. The toric quotient by \(G_{24}\) has lattice width 3. | **In flight:** upstream branch `stack/pr221-cycle-q2` adds a `solution.md`. Coordinate before touching it. |
| `four-plus-two-residual` | #120 | The two-reflection theorem (Ch. 3B) supports all invariant \(H^4\). | Any curve correspondence spanning the anti-invariant \(H^3\) of the weighted threefold \(u^3=H\) (\(h^{2,1}=48\)) needs a parameter curve of genus ≥ 48. The core keeps only one diagonal involution. | An exact suspension reduction \(2A-2dH(v)+3H\operatorname{div}v=0\) isolates the problem on the 4-variable core. |
| `quadric-contact-residual` | #125 | §3D.8 (paired exchange) proves existence non-constructively: it uses the Hodge conjecture on a second sextic \(X'\) plus density. | The paired-cubic Jacobian surface has χ(O)=3 and Kodaira dimension 1, so it is not uniruled. The conic-bundle base is a quintic surface, or a three-nodal quartic K3 on a slice. | A conic-bundle structure in the invariants \((u_i,p_i)\). Compare §3B.11. |

### B. Residual sectors

The same or a closely related variety already has a solved sector, and the target
carries the "wrong" character.

| Target | Discussion | Symmetry and sector | What is missing | Leads |
| --- | --- | --- | --- | --- |
| `cyclic-triple-residual` | #155 | \(G_{18}\cong C_6\times C_3\); the eigenspace has type (0,41,119,26,0). | Cubic pullback supplies only sources with \(X_2\)-exponent ≡ 2 (mod 3); the target has exponent 0. Dividing by \(X_2^2\) introduces poles and a covering-field character. | The free-line Fano quotient and its incidence are already explicit. **Hypothesis:** test the target on the lifted curves with the §5.5 finite exactness test (§6.5). |
| `reflection-coupled-residual` | #154 | Under \(\tau\) and \(Y_2\mapsto-Y_2\), the target is in the \(Y_2\)-even sector, of type (126,492,126). | The quotient \(u=Y_2^2\) has density \(2Y_2\), so it serves only odd targets. | The weighted quartic, deformation system, and quotient primitive are explicit. **Hypothesis:** branch-locus correspondences, or the §5.5 test on double-cover curves. |
| `one-reflection-skew` | #127 | \(G_{24}\); character \((1,1,\zeta_6)\); eigenspace of type (0,44,77,3,0). | The fibre threefolds' 3-forms carry the target's character, which puts it outside Ch. 3D's criteria. The hidden exchanges are broken. | Replacement covers. The lifting and transport machinery only helps once a serving geometry exists. |
| `off-paired` | #57 | \(G_{36}\) fixes the residue. The sector has type (0,22,50,22,0), and all six reflections act by −1 on it. | The paired-sign cubic chart fails, and the Fermat cover and cyclic-six theorem do not apply. | The quotient by an order-6 subgroup gives a weighted CY threefold pencil whose 3-form is fixed by a residual order-3 symmetry. **Hypothesis:** compute the §3D.6 sector count (\(n-c>0\)?) for this pencil. |
| `klein-coupled-character` | #116 | An order-7 action with weights (1,2,4,0,0,0); the target has character \(\zeta^3\); eigenspace of type (0,26,123,37,0). | No target-fixing junior element exists; averaging kills the target. | Character-graded searches over \(\mathbb Q(\zeta_7)\) followed by a trace. Correspondences between characters. |
| `diagonal-permutation` | #177 | The target lies in the 2-dim standard representation of \(S_3\), acting by simultaneous pair permutations. | Scalar character transfer cannot reach a 2-dim representation. | Matrix-valued (equivariant) correspondences; transposition and 3-cycle fixed loci as branch strata. |

### C. Structured equations with no theorem yet

New geometry is needed.

| Target | Discussion | Structure | Notes |
| --- | --- | --- | --- |
| `polarized-k3-remainder` | #115 | A block swap of the K3 double plane \(w^2=f\). | The remainder is \(\Lambda^2H^2(K3)_{\rm prim}\), of type (19,172,19), which is Kuga–Satake territory. Probably the hardest open target. |
| `rank-one-biquartic` | #118 | \(f(X)+g(Y)+t\,L(X)^2M(Y)\); a binary-cubic fibration with discriminant \(-d(4b^3+27a^2d)\). | A deformation of the split case. Analyze \(L=0\), \(M=0\), and the discriminant divisors. |
| `short-power-presentation` | #174 | Eight-power presentation: a codimension-2 linear section of the Fermat sixfold in \(\mathbb P^7\). | Diagonal covers and vanishing cycles of the linear sections. |
| `quadratic-cubic-pullback` | #175 | \(F=\sum Q_i^3\): the pullback of the Fermat cubic fourfold. | The target lies outside the determinant-density direction. The pair-cover transfer leads to Kodaira dimension ≥ 1. |
| `eight-plane-configuration` | #176 | Eight explicit planes; \(F\in(Q_0,Q_1,Q_2)\). | The residual surface fibration has \(p_g=19\), so the \(p_g=0\) criterion fails. The plane classes are of type (2,2), not (3,1). |
| `dense-rational` | #59 | All 236 allowed monomials, near Fermat. | A benchmark with no structural promise. Deformation from Fermat needs control of the passage from a formal series to an exact field. |

**If the goal is a solve soon,** group B targets whose quotient geometry is already
explicit look best: `cyclic-triple-residual` and `reflection-coupled-residual`. Only the
character step is missing, and it can be tested with finite computations (§6.5).

**If the goal is the highest value,** the group A targets turn existing book theorems
into formal proofs, but each needs a new constructive mechanism.

**Avoid for now:** `cycle-q2` (work is in flight), `polarized-k3-remainder`, and
`dense-rational`.

---

## 4. Workflow

Each phase ends in a short written artifact in a per-target dossier (§5.2).

### Phase 0: Due diligence (½ day)

- [ ] Read `data/<id>/why.md`, `instance.rec`, and `provenance.json`, plus the linked
      book sections and the solved companion instances.
- [ ] Read the discussion thread, including replies and every page of comments.
- [ ] Look for open PRs against `hodge-computer/sextics` and upstream branches that touch
      `data/<id>/` (`git log upstream --all -- data/<id>`).
- [ ] Note which book theorems claim support for the target's sector, and whether their
      proofs are constructive.

### Phase 1: Symmetry and character ledger (½ day)

- [ ] Find the projective diagonal symmetry group of \(F\) from the exponent lattice
      (Smith normal form of the exponent differences).
- [ ] Find the linear symmetries beyond diagonal ones: Hessian symmetrizer, apolar base
      locus, permutations.
- [ ] For every symmetry, record the characters of \(A\), of the residue volume, and of
      the residue (their product).
- [ ] Compute the Hodge numbers of the target's eigenspace from the character-graded
      Jacobian ring. The `why.md` files now quote many of these; recompute them as a
      check.
- [ ] List the accepted constructions on this variety and its companions, and which
      sector each one serves.
- **Output:** a table of symmetry, character, and served-or-not.

### Phase 2: Geometry triage (1 day)

- [ ] List candidate reductions: paired-sign cubic charts, conic bundles, cyclic and
      Fermat covers, quotients, linear sections, and pencils of threefolds.
- [ ] For each fibration, compute the fibre invariants χ(O), \(p_g\), and κ, restricted
      to the target's character (extend `scripts/audit-jacobian-scope.py`).
- [ ] Decide:

| Fibres in the target sector | Route |
| --- | --- |
| Rational or uniruled | Integrate along rational curves, using `PolynomialCurves` lifting and a trace. |
| \(p_g=0\), κ ≥ 0 (Dolgachev or elliptic type) | Pass to the Jacobian or relative Jacobian; transport with `ResidueTransport`. |
| Weighted CY threefold pencil | Hodge-locus sector counts (Ch. 3D, §3D.6). Existence only, unless curves can be made explicit. |
| \(p_g>0\) in the sector | This fibration cannot serve the target; look for another. |

### Phase 3: Galois screen (an hour or two)

For every finite configuration the route depends on (sections, lines, conics, the seeds
of an incidence system):

- [ ] Specialize the passive parameters to random rationals and solve over Q with
      `msolve -P 1`.
- [ ] Factor the eliminating polynomial with python-flint.
- [ ] Collect Frobenius cycle types modulo 50–100 primes.
- **Read the result as follows:**
  - An irreducible eliminant at one good specialization proves a single generic Galois
    orbit, because specializing can only split a factorization further.
  - Cycle types identify the group. In `cycle-three`, the 240 sections showed \(30^8\)
    (a Coxeter element), \(24^{10}\), and \(7^{34}1^2\), which points to \(W(E_8)\).
  - **Large orbit:** no explicit support or curve family exists over a small field, so
    go non-constructive (one nonsingular seed, étale lifting, normalized trace).
  - **Small orbit:** explicit sections or curves exist, so try the prescribed-support
    linear algebra (Phase 5).
- **Don't** try to solve an incidence system for all its solutions. Count them first:
  degree-1 del Pezzo fibres have 2,160 conic classes.

### Phase 4: Choose the route

```text
Is the target's sector served by a known theorem?
├─ yes, constructively   → reuse its construction (only targets and character change)
├─ yes, Hodge-theoretic  → group A: find an algebraic realization of the same support
└─ no                    → serving geometry first (Phase 2), then:
     Galois orbit small? ─ yes → explicit support → linear algebra for (D, N, Q)
                          └ no  → seed + lifting + trace (+ transport corrections)
```

### Phase 5: Discovery computations

- Use exact or modular linear algebra for fixed-denominator systems (book §5.1–5.3). At
  q=1 the reduced support \(D_{\rm red}\) suffices (prescribed-support theorem), so only
  the support has to be guessed.
- Use Jacobian-ring kernels of multiplication by \(A\) (the first normalized equation of
  §5.2) as a quick sieve on denominator degrees. This is `rbajaj5`'s method.
- Use Griffiths–Dwork and Gauss–Manin reductions for fibrewise sources. A rational
  solution of the connection means the support is only the fibres. Inconsistency means
  sections must enter, so go back to Phase 3.
- Find seeds with msolve modulo a small prime (under a second). Check the Jacobian
  determinant, and choose the seed that makes the Lean side cheapest (small entries,
  sparse inverse).
- Keep coefficients generic (the family's own parameters) and specialize through ring
  homomorphisms at the end.

### Phase 6: Formalization design

Do this before writing the full proof.

- [ ] List every large polynomial identity the proof needs. Estimate its term counts
      with python-flint.
- [ ] For each identity, pick a checking strategy (§6.2–6.3) and prototype the largest
      one first. Measure the compile time, peak memory, and `leanchecker` replay time.
- [ ] Prefer shared modules; put new reusable machinery in a separate maintenance PR
      first.
- [ ] Budget against 14 GB per process, remembering that importing Mathlib alone costs
      about 3.1 GB.

### Phase 7: Verification and submission

- [ ] Prune modules not reachable from `Solution.lean`: CI compiles every submitted file.
- [ ] Run the local mirror of the protected verifier
      ([`verify_local.py`](cycle-three-compaction/scripts/verify_local.py)) under the
      pinned toolchain. Record the time and memory of each module.
- [ ] Run `scripts/proposals/check-explanation.py` on `solution.md`.
- [ ] Follow the CONTRIBUTING.md file boundary and the PR template. Credit prior results
      and collaborators.

### Phase 8: Share

- [ ] Post or update one discussion comment per approach: the result, its scope, the
      search space, and a link to the scripts. Negative results must state their bounds.

---

## 5. New ways of working

### 5.1 A screening sprint across all open targets

Run Phases 1–3 on all 15 open targets before committing to one, using one script per
screen:

- the character ledger;
- the fibration invariants;
- the Galois orbits;
- the Jacobian-ring kernel dimensions.

Rank the targets by two scores: whether a serving geometry exists, and whether the
heaviest Lean identity looks cheap. This costs about a day and replaces picking a target
on intuition.

### 5.2 A dossier per target

Keep `research/<target>/README.md` up to date, holding:

- the ledger, the triage table, and the screen outputs;
- hypotheses with their status (open, refuted, confirmed);
- links to scripts and logs.

New agents and collaborators can then start from the dossier instead of the chat history.

### 5.3 Hypothesis-driven experiments

Write each idea as a falsifiable test with a stopping rule. For example: "The lifted
free lines on the `cyclic-triple-residual` Fano quotient have genus g, and the target
density is exact on them (by the §5.5 test)." Record a refutation as carefully as a
success; it bounds the search for everyone.

### 5.4 Parallel exploration with adversarial checks

Independent screens and candidate mechanisms parallelize well across subagents or a
workflow. Give each a narrow brief and the exact files to read. Before believing any
positive claim, have a second agent check it adversarially, reading the code rather than
the conclusion. Keep the fan-out small: every agent costs tokens.

### 5.5 Treat certificates as engineering artifacts

- Generators live in the repo and reproduce the Lean data byte for byte.
- Each certificate records its gradings, digit bounds, and coefficient masses, all
  validated outside Lean first.
- Track per-module compile, replay, and peak-memory numbers in the dossier from day one.

### 5.6 Tooling up first when it pays across targets

The graded Kronecker certificate (§6.2) and a library of screen scripts benefit every
target. Land them as small maintenance PRs before the next solve, so the solve itself
stays short.

### 5.7 Make environments reproducible

The container is rebuilt between sessions, and Singular, msolve, python-flint, elan, and
the Mathlib cache all disappear. Save them as a setup script on the cloud environment
(see §7) so every session starts ready.

---

## 6. Technique cards

### 6.1 Galois-orbit screen

```sh
# from research/cycle-three-compaction/cas (needs msolve and python-flint)
python3 gen_sections.py <c1> <c2> 0 sys.ms     # write the system over Q at a specialization
msolve -f sys.ms -o sys.out -P 1 -t 4          # rational parametrization
python3 factor_elim.py sys.out                 # factor degrees = Galois orbit sizes
python3 frobenius.py <c1> <c2>                 # cycle types modulo 60 primes
```

Adapt `gen_sections.py` to the configuration at hand (sections of a Weierstrass model,
lines on a cubic surface, conics, seeds of an incidence system).

### 6.2 Graded whole-identity Kronecker certificate

From [`NineKronecker.lean`](cycle-three-compaction/lean/NineKronecker.lean) and
[`FullSyzygy.lean`](cycle-three-compaction/lean/FullSyzygy.lean):

1. Find a nonnegative \(\mathbb N^k\) weight grading under which the identity is
   homogeneous. Keeping the homogeneous coordinate \(z\) helps; for `cycle-three` it gave
   four components.
2. The variables whose exponents the grade determines are set to 1. The others are
   Kronecker-encoded, with strides given by their degree bounds (radix 19 there). In
   `cycle-three` this left 5 of 9 variables encoded.
3. Choose \(B=2^k\) above the coefficient mass, computed from `splitEval` at all ones
   (\(2^{50}\) there).
4. Write the polynomials in homogeneous Horner form, then run one `decide +kernel` on
   both sides as single big integers. The kernel uses GMP, so this is only a few dozen
   large multiplications.
5. Result for `cycle-three`: 1.65 MB to 26 KB, about 21 CPU-minutes to 40 s, and peak
   memory under the Mathlib baseline plus 2 GB.

The shared `KroneckerCheck` instead gets injectivity from per-variable degree bounds and
a grid sweep. Grading removes variables altogether, so the two methods combine well.

### 6.3 Replacing large `ring` or `simp` proofs

Instantiate the certificate's soundness theorem (it holds in every commutative ring) at
the ring where the statement lives. For example, use `MvPolynomial (Fin 2) K` with
\(x,y,z\mapsto X_0,X_1,1\) and coefficients mapped to `C a`. Prove the derivative or
evaluation bridges with small `ring` calls. In `cycle-three`, `CubicDifferential` went
from 324 s and 8.9 GB to 110 s and 5.1 GB.

### 6.4 Lean kernel pitfalls

- The kernel keeps every intermediate result for the whole declaration. Symbolic product
  trees or one giant `decide` reached 12 GB or more. Split the work, or compute on
  integers rather than syntax trees.
- Concrete huge exponents inside proofs (`X^130321`) make the kernel unfold `npow` until
  it hits deep recursion. Keep strides symbolic (`X^(strides i)`).
- One `omega` call over many weight equations ran out of memory. Use small helper lemmas.
- `leanchecker` replay re-runs every `decide`; for the original `cycle-three` it cost as
  much as compiling.
- CI compiles unused files, so prune them.

### 6.5 Testing rational sources on a known quotient

**This is a hypothesis to test, aimed at group B.** Many residual targets fail because
the transported source acquires poles, or a covering-field character, on an otherwise
serving quotient. Proposed test:

1. Lift the explicit curve family of the solved sector to the cover on which the target's
   source is a function.
2. Compute the genus of the lifted curves, and the target's differential along them.
3. Apply the §5.5 exactness test: \(\omega\in d(L(D_\omega))\), which is finite linear
   algebra over the passive field. The `QuadraticCovers` residue machinery and §3B.11's
   global residue cancellation show how to control spurious poles.
4. **If the differential is exact:** a primitive exists along the family. Assemble it with
   the existing transport and trace.
5. **If it is not:** the regular differentials of the lifted curves are the obstruction,
   which is itself a precise negative result worth posting.

### 6.6 Prescribed-support linear algebra

At q=1, with a guessed reduced support \(D_{\rm red}\), solve \(N(F)=A\,D\) and
\(D\operatorname{div}N-N(D)\in(F)\) as exact or modular linear systems. Restricting to
the right character divides the number of unknowns by roughly the order of the symmetry
group. Use it only after the Galois screen says a small support is plausible.

---

## 7. Tool setup

Run these in the container. The last two lines assume a copy of `formal/lean`.

```sh
apt-get install -y singular msolve        # macaulay2 is also available
pip install python-flint sympy
curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none
export PATH="$HOME/.elan/bin:$PATH"
(cd formal/lean && lake exe cache get)    # about 7.6 GB; fetches the pinned toolchain and Mathlib
```

Save these as the environment's setup script so they persist across sessions.

---

## 8. Lessons from `cycle-three` (PR #165)

- The mathematics was natural and close to minimal:
  - a paired-sign cubic fibration;
  - a covariant map to the Jacobian;
  - rational elliptic surfaces, which are degree-1 del Pezzo surfaces after blowing down;
  - polynomial curves found from a mod-5 seed, lifted;
  - a trace, then transport with moving-map and Piola corrections.
- 91% of the 2.14 MB proof was one identity checked inefficiently. A variant with graded
  Kronecker certificates is 507 KB. Under the same verifier steps it compiles in 917 s
  instead of 2,188 s, replays in 443 s instead of 1,506 s, and peaks at 5.1 GB instead
  of 8.9 GB.
- CAS tools did not find a better certificate; they explained why none is easy to find.
  The 240 sections form one \(W(E_8)\)-like Galois orbit over the parameter field, so the
  non-constructive lifting is the economical design. Reducing the covariants modulo the
  cubic did not shrink them.
- The original surfaces are Dolgachev-type, so the Jacobian detour is essential: they
  carry no rational-curve families. The `cycle-q2` and `quadric-contact-residual` notes
  meet the same obstruction (a non-uniruled elliptic surface with χ(O)=3).

Details: [cycle-three-compaction/README.md](cycle-three-compaction/README.md).
