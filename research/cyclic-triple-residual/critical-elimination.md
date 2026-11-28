# Direct critical elimination: measured limits

Research checkpoint, 2026-11-28, by OpenAI Codex for
`hodgepodgeconjecture`. This records actual attempts to compute the
complete ordinary critical algebra, not another local-point test.

The fifteen-equation multiplier system was specialized to K=1,L=2.
Its equation term counts are 4,9,16,21,22,15,8,20,20,42,28,17,17,17,17.
The ideal was unsaturated: even a completed basis would still need removal
of the excluded components before yielding the trace.

- SymPy's F5B implementation, grevlex order over F13, did not finish
  within 60 seconds under a 3 GiB address-space limit.
- FLINT 0.9.0 was installed in a local isolated dependency directory.
  Its integer-polynomial Buchberger implementation returned an incomplete
  generating set of length 57 when given limits of 200 basis elements,
  2,000 terms per polynomial, and 512 coefficient bits. The internal
  growth limit was reached in approximately 0.028 seconds.
- With limits increased to 1,000 basis elements, 50,000 polynomial terms,
  and 4,096 coefficient bits, that implementation did not finish before
  the external 60-second timeout. The process exited with code 124.

These tests use different arithmetic and are not a speed comparison.
Neither is an F4 implementation. Neither yielded a complete Groebner
basis, a quotient dimension, or the global trace. No conclusion about
existence of a primitive follows from the failures.

[critical-elimination-attempts.json](critical-elimination-attempts.json)
records the observed outcomes. A reusable bounded runner is supplied:

```sh
python3 research/cyclic-triple-residual/critical_elimination_probe.py --backend sympy --seconds 60
# Requires python-flint==0.9.0 in the active Python environment:
python3 research/cyclic-triple-residual/critical_elimination_probe.py --backend flint --seconds 60
```

The parent process imposes the timeout even while compiled code runs.
Its one-second timeout handling was smoke-tested. The recorded 60-second
runs used the preceding prototype implementations of these same systems
and algorithm settings; they are not represented as completed runs of
the packaged runner.

The next implementation to test should exploit the already established
three-slope elimination, or use a sparse modular F4 solver. Simply raising
the current limits has not been justified. No AWS resources were used.
The mathematical obligations concerning saturation, omitted boundary
terms, and the holomorphic remainder remain unchanged.

See [principal-parts.md](principal-parts.md) for the precise trace and
[README.md](README.md#fixed-inputs-and-prior-work) for attribution.
No protected files, discussion posts, or PRs were changed.
