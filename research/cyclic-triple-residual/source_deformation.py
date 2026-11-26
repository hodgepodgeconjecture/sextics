"""Audit the exact source-deformation current on the fixed incidence component.

F_s=F-s*A is an auxiliary differentiation device; the target remains s=0.
The determinant currents compare its residue derivative with the original
normal representative. The proof is in README; finite-precision checks are
independent implementation checks, not generic vanishing proofs or Lean.
"""
import argparse
import json
from pathlib import Path
import sympy as sp
from line_exactness import (compile_poly, derivative, evaluate, inverse, mv,
                            add, mul, scale, power, diff, det_poly)
from source_root_chart import Adic, inv_adic, build


def audit(p, n):
    Adic.p = p
    Adic.cap = 2*n
    modulus = p**n
    Q, D, DD = build('X1')
    ri = next(i for i in range(p) if (i*i+1) % p == 0)
    point = [ri, 0, 0, ri, 0, 0, 0, p, 1, 2*p, 3*p]
    for _ in range(n.bit_length()+2):
        jac = [[evaluate(d, point, modulus) for d in row[:7]] for row in D]
        delta = mv(inverse(jac, modulus), [evaluate(q, point, modulus) for q in Q], modulus)
        point[:7] = [(x-v) % modulus for x, v in zip(point[:7], delta)]
    assert all(evaluate(q, point, modulus) == 0 for q in Q)
    point = [Adic(v, n) for v in point]
    jac = [[Adic(evaluate(d, point, 1)) for d in row] for row in D]
    ji = inv_adic([row[:7] for row in jac])
    hess = [[[Adic(evaluate(dd, point, 1)) for dd in row] for row in qq] for qq in DD]
    directions = [mv(ji, [7, 0, 0, 0, 0, 0, 0], 1)+[Adic(0)]*4]
    # X1 is anchored, hence the k1 current is identically zero.
    for j in [7, 9, 10]:
        directions.append(mv(ji, [-row[j] for row in jac], 1)
                          +[Adic(int(j == l)) for l in range(7, 11)])
    def second(u, v):
        forcing = [sum(h[i][j]*u[i]*v[j] for i in range(11) for j in range(11)) for h in hess]
        return mv(ji, [-x for x in forcing], 1)+[Adic(0)]*4
    def curve(v):
        return [[v[0], v[1]], [v[7], v[2]], [v[9], v[3]], [v[10], v[4], v[5], v[6]]]
    def matrix(cols):
        return [[col[i] for col in cols] for i in range(4)]
    def replace(m, j, col):
        return [[col[i] if l == j else m[i][l] for l in range(4)] for i in range(4)]
    def determinant_derivative(m, dm):
        out = [Adic(0)]
        for j in range(4):
            out = add(out, det_poly(replace(m, j, [dm[i][j] for i in range(4)]), 1), 1)
        return out
    def zero(polys):
        values = [Adic(v) for poly in polys for v in poly]
        assert all(v == 0 for v in values)
        return min(v.n for v in values)

    xi, a, b, c, d, e, f, k0, k1, k2, k3 = point
    t = sp.Symbol('t')
    vs = sp.symbols('xi a b c d e f k0 k1 k2 k3')
    Apoly = sp.Poly((vs[7]+vs[2]*t)**3*(vs[0]+vs[1]*t)**3, t)
    Aq = [compile_poly(Apoly.nth(i), vs) for i in range(7)]
    AD = [[derivative(q, j) for j in range(11)] for q in Aq]
    A = [Adic(evaluate(q, point, 1)) for q in Aq]
    U = mv(ji, [7*v for v in A], 1)+[Adic(0)]*4
    # Mixed source-deformation derivatives: Q_y U=7*A, including A's jet.
    Ud = []
    Aderivatives = []
    for v in directions:
        Av = [sum(Adic(evaluate(dq, point, 1))*vj for dq, vj in zip(row, v)) for row in AD]
        Aderivatives.append(Av)
        correction = mv(ji, [7*x for x in Av], 1)+[Adic(0)]*4
        Ud.append([x+y for x, y in zip(second(U, v), correction)])

    M = matrix([curve(v) for v in directions])
    J = det_poly(M, 1)
    Jh = determinant_derivative(M, matrix([curve(second(v, directions[0])) for v in directions]))
    Js = determinant_derivative(M, matrix([curve(v) for v in Ud]))
    expected_J = 1/(6*xi**5)
    expected_Jh = -5/(36*xi**11)
    expected_Js = -k0**3/(18*xi**8)
    determinant_precision = zero([add(J, [-expected_J], 1),
                                  add(Jh, [-expected_Jh], 1),
                                  add(Js, [-expected_Js], 1)])
    UC = curve(U)
    S = [det_poly(replace(M, j, UC), 1) for j in range(1, 4)]
    Sd = []
    for j in range(1, 4):
        Mj = matrix([curve(second(v, directions[j])) for v in directions])
        Sd.append(determinant_derivative(replace(M, j, UC), replace(Mj, j, curve(Ud[j]))))

    # Cramer's rule gives U=A*V+sum_i (S_i/J)*K_i in physical coordinates.
    cramers = []
    for row in range(4):
        rhs = mul(mul(A, J, 1), M[row][0], 1)
        for j in range(3): rhs = add(rhs, mul(S[j], M[row][j+1], 1), 1)
        cramers.append(add(mul(J, UC[row], 1), scale(rhs, -1, 1), 1))
    cramer_precision = zero(cramers)
    AJh = add(mul(Aderivatives[0], J, 1), mul(A, Jh, 1), 1)
    divergence = AJh
    for v in Sd: divergence = add(divergence, v, 1)
    divergence_precision = zero([add(Js, scale(divergence, -1, 1), 1)])

    Z = [k3, d, e, f]
    Zh = curve(directions[0])[3]
    Zs = UC[3]
    oldR = add(scale(mul(Z, AJh, 1), 3, 1), scale(mul(mul(A, J, 1), Zh, 1), -2, 1), 1)
    newR = add(scale(mul(Z, Js, 1), 3, 1), scale(mul(J, Zs, 1), -2, 1), 1)
    currentR = [Adic(0)]
    for j in range(3):
        term = add(scale(mul(Z, Sd[j], 1), 3, 1),
                   scale(mul(S[j], curve(directions[j+1])[3], 1), -2, 1), 1)
        currentR = add(currentR, term, 1)
    transport_precision = zero([add(newR, scale(add(oldR, currentR, 1), -1, 1), 1)])

    def L(v):
        return add(scale(mul(Z, diff(v), 1), 3, 1), scale(mul(diff(Z), v, 1), -2, 1), 1)
    cols = [L([Adic(0)]*j+[Adic(1)]) for j in range(3)]+[Z, [Adic(0)]+Z]
    mat = [[col[i] if i < len(col) else Adic(0) for col in cols] for i in range(5)]
    sol = mv(inv_adic(mat), newR+[Adic(0)]*(5-len(newR)), 1)
    normal_precision = zero([add(newR, scale(add(L(sol[:3]), mul(Z, sol[3:], 1), 1), -1, 1), 1)])
    aa, bb, cc, dd = sp.symbols('aa bb cc dd')
    disc = sp.discriminant(aa*t**3+bb*t*t+cc*t+dd, t)
    dc = compile_poly(disc, [aa, bb, cc, dd])
    zreverse = list(reversed(Z))
    Delta = Adic(evaluate(dc, zreverse, 1))
    Delta_s = sum(Adic(evaluate(derivative(dc, j), zreverse, 1))*list(reversed(Zs))[j] for j in range(4))
    predicted_C0 = 3*expected_Js-expected_J*Delta_s/(2*Delta)
    scalar_precision = zero([[sol[3]-predicted_C0], [sol[4]]])
    assert sol[3] != 0  # This sample retains the compact elliptic obstruction.
    def record(v):
        return {'valuation':v.val(), 'absolute_precision':v.n,
                'normalized_residue_mod_p':int(v.q/p**v.val()) % p if v.q and v.val() >= 0 else None,
                'zero_to_precision':v == 0}
    return {'prime':p, 'input_precision':n, 'target_source':'X0^3*Y0^3',
            'auxiliary_deformation':'F-s*A; all reported source identities at s=0',
            'current_parameters':['k0','k2','k3'], 'current_numerator_degrees':[len(v)-1 for v in S],
            'old_R_degree':len(oldR)-1, 'new_R_degree':len(newR)-1,
            'checked_to':{'constant_determinants':determinant_precision, 'cramers_rule':cramer_precision,
                          'volume_divergence':divergence_precision, 'full_source_transport':transport_precision,
                          'normal_form_reconstruction':normal_precision, 'scalar_discriminant_formula':scalar_precision},
            'new_compact_remainder':record(sol[3]), 'new_residue_remainder':record(sol[4]),
            'scope':'Exact generic current formulas are proved in README; these adic samples independently check the implementation. No target certificate or Lean verification.'}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--prime', type=int, default=101)
    ap.add_argument('--precision', type=int, default=30)
    ap.add_argument('--output', type=Path, default=Path(__file__).with_name('source-deformation-101.json'))
    args = ap.parse_args()
    result = audit(args.prime, args.precision)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
