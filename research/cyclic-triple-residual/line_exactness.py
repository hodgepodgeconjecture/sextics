"""Exact modular/adic audit of the fixed cyclic-triple-residual weighted-line route.

Uses only Python and SymPy. No numerical roots or floating-point arithmetic.
The local incidence is the accepted cyclic-triple-fano family's anchored chart.
This is research, not a Lean certificate or a nonexistence proof for the target.
"""
import argparse
import itertools
import json
import time
from pathlib import Path
import sympy as sp


def inverse(a, modulus):
    n = len(a)
    m = [[int(v) % modulus for v in row] + [int(i == j) for j in range(n)]
         for i, row in enumerate(a)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if sp.gcd(m[i][j], modulus) == 1)
        m[j], m[pivot] = m[pivot], m[j]
        v = pow(m[j][j], -1, modulus)
        m[j] = [v*x % modulus for x in m[j]]
        for i in range(n):
            if i != j:
                v = m[i][j]
                m[i] = [(x-v*y) % modulus for x, y in zip(m[i], m[j])]
    return [row[n:] for row in m]


def mv(a, x, modulus):
    return [sum(u*v for u, v in zip(row, x)) % modulus for row in a]


def compile_poly(expr, variables):
    return [(int(c), tuple(e)) for e, c in sp.Poly(expr, *variables).terms()]


def evaluate(terms, point, modulus):
    powers = [[pow(x, n, modulus) for n in range(7)] for x in point]
    ans = 0
    for c, exps in terms:
        for j, e in enumerate(exps):
            if e:
                c = c*powers[j][e] % modulus
        ans += c
    return ans % modulus


def derivative(terms, j):
    ans = []
    for c, e in terms:
        if e[j]:
            ee = list(e)
            ee[j] -= 1
            ans.append((c*e[j], tuple(ee)))
    return ans


def add(a, b, mod):
    out = [0]*max(len(a), len(b))
    for i, v in enumerate(a): out[i] += v
    for i, v in enumerate(b): out[i] += v
    return trim([v % mod for v in out])


def scale(a, c, mod): return trim([x*c % mod for x in a])


def trim(a):
    while len(a) > 1 and a[-1] == 0: a.pop()
    return a


def mul(a, b, mod):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i+j] += x*y
    return trim([x % mod for x in out])


def power(a, n, mod):
    out = [1]
    for _ in range(n): out = mul(out, a, mod)
    return out


def diff(a): return [i*a[i] for i in range(1, len(a))] or [0]


def det_poly(m, mod):
    n = len(m)
    out = [0]
    for perm in itertools.permutations(range(n)):
        sign = (-1)**sum(perm[i] > perm[j] for i in range(n) for j in range(i+1, n))
        term = [sign]
        for i in range(n): term = mul(term, m[i][perm[i]], mod)
        out = add(out, term, mod)
    return out


def rank(a, p):
    a = [[int(v) % p for v in row] for row in a]
    row = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(row, len(a)) if a[i][col]), None)
        if pivot is None: continue
        a[row], a[pivot] = a[pivot], a[row]
        v = pow(a[row][col], -1, p)
        a[row] = [x*v % p for x in a[row]]
        for i in range(row+1, len(a)):
            v = a[i][col]
            if v: a[i] = [(x-v*y) % p for x,y in zip(a[i], a[row])]
        row += 1
        if row == len(a): break
    return row


def main(p, precision, output):
    start = time.monotonic()
    modulus = p**precision
    t = sp.Symbol('t')
    xi, a, b, c, d, e, f, k0, k1, k2, k3 = variables = sp.symbols('xi a b c d e f k0 k1 k2 k3')
    y0, y1, x0, x1, y2, z = xi+a*t, t, k0+b*t, k1, k2+c*t, k3+d*t+e*t*t+f*t**3
    H7 = 7*(y0**6+y1**6+x0**6+x1**6+y2**6+z*z) + 5*(
        x0**3*x1*y0*y1+x0**3*z+x0*x0*y0*y1*y2*y2+
        x0*x1*y1**3*y2+x0*z*y1*y2+x1*z*y0*y1)
    poly = sp.Poly(sp.expand(H7), t)
    Q = [compile_poly(poly.nth(i), variables) for i in range(7)]
    D = [[derivative(q, j) for j in range(11)] for q in Q]
    DD = [[[derivative(row[j], l) for l in range(11)] for j in range(11)] for row in D]
    root_i = next(i for i in range(p) if (i*i+1) % p == 0)
    point = [root_i, 0, 0, root_i, 0, 0, 0, p, 1, 2*p, 3*p]
    for iteration in range(precision.bit_length()+2):
        values = [evaluate(q, point, modulus) for q in Q]
        jac = [[evaluate(dq, point, modulus) for dq in row[:7]] for row in D]
        ji = inverse(jac, modulus)
        correction = mv(ji, values, modulus)
        point[:7] = [(x-v) % modulus for x,v in zip(point[:7], correction)]
    assert all(evaluate(q, point, modulus) == 0 for q in Q)
    jac_all = [[evaluate(dq, point, modulus) for dq in row] for row in D]
    ji = inverse([row[:7] for row in jac_all], modulus)
    eps = mv(ji, [7,0,0,0,0,0,0], modulus)+[0]*4
    directions = [eps]
    for j in [7,9,10]:
        dy = mv(ji, [-row[j] for row in jac_all], modulus)
        directions.append(dy+[int(j==l) for l in range(7,11)])
    hess = [[[evaluate(ddq, point, modulus) for ddq in row] for row in qq] for qq in DD]
    second = []
    for direction in directions:
        forcing = [sum(h[i][j]*eps[i]*direction[j] for i in range(11) for j in range(11)) % modulus for h in hess]
        ddy = mv(ji, [-v for v in forcing], modulus)
        second.append(ddy+[0]*4)
        # Differentiate the actual seven coefficient equations, not just the curve value.
        assert all((v+u) % modulus == 0 for v,u in zip(mv([row[:7] for row in jac_all],ddy,modulus),forcing))
    def curve_dir(v):
        return [[v[0],v[1]], [v[7],v[2]], [v[9],v[3]], [v[10],v[4],v[5],v[6]]]
    columns = [curve_dir(v) for v in directions]
    m = [[col[i] for col in columns] for i in range(4)]
    J = det_poly(m, modulus)
    Jdot = [0]
    for j in range(4):
        mm = [[v[:] for v in row] for row in m]
        cc = curve_dir(second[j])
        for i in range(4): mm[i][j] = cc[i]
        Jdot = add(Jdot, det_poly(mm,modulus),modulus)
    # Degree argument in README: the evaluation determinant is time-constant.
    # Retain the independently assembled determinant as a check of this shortcut.
    expected_J=pow(6*pow(point[0],5,modulus),-1,modulus)
    expected_Jdot=-5*pow(36*pow(point[0],11,modulus),-1,modulus) % modulus
    assert J==[expected_J] and Jdot==[expected_Jdot]
    X = [point[7],point[2]]; Y = [point[0],point[1]]
    Z = [point[10],point[4],point[5],point[6]]
    Zdot = curve_dir(eps)[3]
    A = mul(power(X,3,modulus),power(Y,3,modulus),modulus)
    Adot = scale(add(mul(mul(power(X,2,modulus),power(Y,3,modulus),modulus),curve_dir(eps)[1],modulus),
                      mul(mul(power(X,3,modulus),power(Y,2,modulus),modulus),curve_dir(eps)[0],modulus),modulus),3,modulus)
    AJ = mul(A,J,modulus)
    AJdot = add(mul(Adot,J,modulus),mul(A,Jdot,modulus),modulus)
    R = add(scale(mul(Z,AJdot,modulus),3,modulus),scale(mul(AJ,Zdot,modulus),-2,modulus),modulus)
    def valuation(v):
        if v == 0: return precision
        n=0
        while v % p == 0: n+=1; v//=p
        return n
    vz=min(map(valuation,Z)); vr=min(map(valuation,R))
    Zbar=[v//p**vz % p for v in Z]; Rbar=[v//p**vr % p for v in R]
    assert len(Zbar)==4 and Zbar[-1] != 0
    assert sp.gcd(sp.Poly.from_list(Zbar[::-1],t,modulus=p),sp.Poly.from_list(diff(Zbar)[::-1],t,modulus=p)).degree()==0
    # The constant-determinant lemma sharpens deg R <=9 and deg s <=7.
    degree_bound=7
    row_count=degree_bound+3
    operator=[]
    for n in range(degree_bound+1):
        monomial=[0]*n+[1]
        col=add(scale(mul(Zbar,diff(monomial),p),3,p),scale(mul(diff(Zbar),monomial,p),-2,p),p)
        operator.append(col+[0]*(row_count-len(col)))
    matrix=[list(row) for row in zip(*operator)]
    rhs=Rbar+[0]*(row_count-len(Rbar))
    rk=rank(matrix,p); aug=rank([row+[v] for row,v in zip(matrix,rhs)],p)
    # Reduce the FULL adic source, not just its leading residue-field term.
    # Division by p^vr loses vr known digits; use the smaller common precision.
    reduced_precision=precision-max(vz,vr)
    reduced_modulus=p**reduced_precision
    zn=[v//p**vz % reduced_modulus for v in Z]
    rn=trim([v//p**vr % reduced_modulus for v in R])
    def operator_poly(s):
        return add(scale(mul(zn,diff(s),reduced_modulus),3,reduced_modulus),
                   scale(mul(diff(zn),s,reduced_modulus),-2,reduced_modulus),reduced_modulus)
    reduced=rn[:]; sn=[0]
    while len(reduced)>5:
        n=len(reduced)-3
        coeff=reduced[-1]*pow(3*zn[3]*(n-2),-1,reduced_modulus) % reduced_modulus
        term=[0]*n+[coeff]
        sn=add(sn,term,reduced_modulus)
        reduced=add(reduced,scale(operator_poly(term),-1,reduced_modulus),reduced_modulus)
    cols=[operator_poly([0]*n+[1]) for n in range(3)]+[zn,[0]+zn]
    small_matrix=[[col[i] if i<len(col) else 0 for col in cols] for i in range(5)]
    small_solution=mv(inverse(small_matrix,reduced_modulus),reduced+[0]*(5-len(reduced)),reduced_modulus)
    sn=add(sn,small_solution[:3],reduced_modulus)
    reconstructed=add(operator_poly(sn),mul(zn,small_solution[3:],reduced_modulus),reduced_modulus)
    assert add(reconstructed,scale(rn,-1,reduced_modulus),reduced_modulus)==[0]
    report=dict(prime=p,precision=precision,coupling='5/7',endpoint=[p,1,2*p,3*p],
                incidence_point=point,incidence_equations_mod_pN=True,implicit_second_derivatives_checked=True,
                J_degree=len(J)-1,R_degree=len(R)-1,Z_valuation=vz,R_valuation=vr,
                constant_determinant_shortcut_checked=True,
                normalized_Z_mod_p=Zbar,normalized_R_mod_p=Rbar,
                primitive_degree_bound=degree_bound,operator_shape=[row_count,degree_bound+1],rank=rk,augmented_rank=aug,
                full_normal_form_reconstruction=True,normalized_remainder_precision=reduced_precision,
                normalized_remainders_mod_pN=small_solution[3:],
                normalized_remainders_mod_p=[v % p for v in small_solution[3:]],
                elapsed_seconds=round(time.monotonic()-start,3),
                scope='Local p-adic branch through accepted seed; bounded time-direction equation only. No Lean proof.')
    output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='incidence_point'},indent=2),flush=True)


if __name__ == '__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--prime',type=int,default=101); ap.add_argument('--precision',type=int,default=8)
    ap.add_argument('--output',type=Path,default=Path(__file__).with_name('result.json')); args=ap.parse_args()
    main(args.prime,args.precision,args.output)
