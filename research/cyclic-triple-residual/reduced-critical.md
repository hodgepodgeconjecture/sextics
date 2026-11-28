# Reduced critical equations and the first sparse-solver runs

Research checkpoint, 2026-11-28, by OpenAI Codex for
`hodgepodgeconjecture`. The global trace remains uncomputed.

## Five-variable presentation

Fix K=1,L=2 and work modulo 13. In the endpoint model write D0=w^2.
The three slope equations P4,P5,P6 involve only (a,b,c,xi), since all
powers of w in these equations are even. The reconstructed cubic's
jth coefficient has denominator D0^j. Its discriminant can therefore
be computed in the quadratic algebra with relation w^2=D0, initially
with denominator D0^6. Exact polynomial division cancels one common D0
factor in this specialization, giving

\[
\Delta=(U+wV)/D_0^5.
\]

For the 3-by-4 derivative matrix of P4,P5,P6 in (a,b,c,xi), let N_i be
the signed maximal minors, chosen so N_xi is the three-slope Jacobian.
Then the vector N annihilates all three incidence equations. For an
expression (U+wV)/D0^q, its derivative numerator along this vector is

\[
C_0=\sum_i N_i(D_0U_i-q(D_0)_iU),\qquad
C_1=\sum_i N_i(D_0V_i+(1/2-q)(D_0)_iV).
\]

Thus Delta_xi=(C0+w*C1)/(D0^(q+1)*N_xi) on the ordinary chart.
Here a further common D0^3 factor cancels from C0,C1. Sequential division
by the three incidence equations reduces their size, with every division
identity checked exactly. These are valid representatives modulo the
incidence ideal, not a claim to have computed a Groebner normal form.

The resulting five equations are D0-w^2=0, P4=P5=P6=0, C0+w*C1=0.
They require localization away from D0, N_xi, Delta and the second
Delta derivative to recover exactly the ordinary critical locus.

| Measurement | Result |
| --- | --- |
| Slope polynomial term counts | 211, 437, 854 |
| Slope total degrees (including xi) | 22, 29, 36 |
| Discriminant numerator term counts U,V | 1,039; 964 |
| Critical coefficient terms before cancellation | 93,384; 89,742 |
| After D0 cancellation and incidence reductions | 50,149; 47,528 |
| Final critical coefficient total degrees | 112; 106 |

The reconstructed Delta_xi is checked against both existing witnesses:
it is zero at the ordinary critical point and four at the noncritical
point modulo 13. The derivation gives the algebraic identity; these two
checks verify its implementation. Eliminating variables has produced
substantial degree growth, so the five-variable model is not automatically
preferable to the fifteen-variable multiplier system.

Reproduce with python-flint==0.9.0 installed:

```sh
python3 research/cyclic-triple-residual/reduced_critical.py
# Optional large solver input, not committed:
python3 research/cyclic-triple-residual/reduced_critical.py --export-input /tmp/reduced-critical.ms
```

The compact generator and [report](reduced-critical.json) are saved;
the large expanded input is reproducible and is not included.

## Sparse modular solver

A local msolve 0.11.0 installation was assembled from the JuliaBinaryWrappers
release packages for msolve 0.1100.0+0, FLINT 301.600.1+0, and OpenBLAS32
0.3.34+2. The binary reports msolve 0.11.0. A smoke test on x^2-1 over F13
returned the expected reduced basis. No AWS resources were used.

The fifteen-variable unsaturated system was run with exact sparse linear
algebra, two threads, at most 1,000 pairs per matrix, a 3 GiB address-space
limit, and a 60-second timeout. It reached degree nine and processed a
53,313-by-226,539 sparse matrix, but did not finish. Peak resident memory
was about 959 MiB. This is solver progress, not a completed basis or trace.

The native F4SAT mode rejected characteristic 13 as unsupported. A further
run used the exact integer input reduced modulo 1,000,003, saturating by
xi*(14*k3+5)*Delta. This removes some excluded loci, not all of them.
The large-prime native run exited with signal 11 after about 59 seconds,
with peak resident memory around 1.63 GiB under the 3 GiB address-space
cap. Its last logged kernel had 452,663 rows and 1,397,481 columns. The
cause of the crash is not established. The run reports are saved in [msolve-runs.json](msolve-runs.json).
Different primes and saturation settings prevent treating their timings
as a controlled speed comparison.

A correctness smoke test then failed: native -S returned x^2-1 for
saturation of (x^2-1) by x-1, instead of the required x+1. A zero exit
code was therefore insufficient. We do not rely on native saturation
output. Adding an inverse variable u with u*(x-1)-1 and using ordinary
F4 returned x+1 and u+500002 over F1000003, the correct result.
The runner now uses this Rabinowitsch construction for --saturate,
introducing a sixteenth variable, rather than native -S.

The inverse-variable system was subsequently tested over F13: a
60-second run with two threads and a 3 GiB cap timed out near 963 MiB
resident memory. A 120-second run with four threads and a 6 GiB cap
also timed out, reaching roughly 1.96 GiB resident memory and additional
degree-nine reductions. Neither produced a complete basis. These runs
used the packaged runner. The larger cap was not exhausted by resident
memory, so these data do not establish a need for EC2.

A portable runner exports the fifteen-variable system and enforces limits:

```sh
python3 research/cyclic-triple-residual/msolve_probe.py --msolve /path/to/msolve --prime 1000003 --saturate --seconds 60 --output /tmp/critical-msolve
```

The reports describe the prototype runs; the packaged runner's timeout
path is tested separately. Neither runner computes a multiplication trace
or asserts that a successful solver exit independently certifies a proof.

The solver source and release instructions were read from the
[msolve v0.11.0 repository](https://github.com/algebraic-solving/msolve/tree/v0.11.0).
Credit for msolve belongs to J. Berthomieu, C. Eder, and M. Safey El Din;
their library is described in ISSAC 2021, pp. 51–58. The packaged binaries
are supplied by JuliaBinaryWrappers. The fixed incidence and source work
is credited in [README.md](README.md#fixed-inputs-and-prior-work).

No complete critical algebra, corrected global pairing, rational current,
or Lean target certificate has been obtained. The remaining work includes
full exceptional-locus handling and the holomorphic remainder check.
