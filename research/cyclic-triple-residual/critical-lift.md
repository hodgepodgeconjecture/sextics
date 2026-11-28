# An actual ordinary critical point and its exact lift

Research checkpoint, 2026-11-28, by OpenAI Codex for
`hodgepodgeconjecture`. This validates the [Lagrange critical system](critical-lagrange.md)
on its intended locus, beyond the earlier noncritical test point.

At K=1,L=2 modulo 13, the eight curve coordinates (xi,a,b,c,d,e,f,k3) are
(4,6,2,12,6,9,3,1). The seven multipliers are (3,0,10,1,11,9,3).
They satisfy all seven incidence equations and all eight stationarity
equations. The full 15-by-15 Jacobian has determinant 4 modulo 13.

The discriminant derivative vanishes, while the remaining local data are:

| Quantity | Value modulo 13 |
| --- | ---: |
| Delta | 3 |
| w0 | 12 |
| Delta_xixi | 12 |
| Delta_dot at fixed xi | 5 |
| 2 Delta_dot / (w0^2 Delta Delta_xixi) | 1 |

Thus this is an ordinary simple critical point in the required open chart,
with a nonzero local trace contribution. Newton iteration lifts all fifteen
coordinates modulo 13^32. Substitution into **all fifteen equations** gives
zero modulo 13^32; every denominator used to evaluate the summand is a unit.

The invertible Jacobian also supplies the characteristic-zero conclusion
by Hensel lifting: an isolated smooth algebraic critical point exists, and
its local summand is nonzero. Allowing K,L to vary gives an etale local
critical family. This establishes that the ordinary critical locus and
its nonzero source contribution are not artifacts of a singular reduction.
The recorded finite precision is a reproducibility check, not an assertion
that a finite list of digits itself is a characteristic-zero equality.

Run `python3 research/cyclic-triple-residual/critical_lift.py`.
[critical-lift.json](critical-lift.json) contains the seed, lift, Jacobian
check, and local summand. The script binds the construction to the protected
target through the existing curve-model audit. It retains the actual
source deformation H-sA throughout.

**This is not the global trace.** Other critical points, including points
not rational over the residue field, can cancel this contribution. No
sum over all closed points, global residue-duality obstruction, rational
primitive, or target Lean certificate is established. The next substantive
computation must recover the complete critical algebra or an equivalent
complete trace calculation; more precision at this single point cannot
settle that question.

The incidence and source work used here is credited in
[README.md](README.md#fixed-inputs-and-prior-work). No protected file,
discussion, or PR is changed by this checkpoint.
