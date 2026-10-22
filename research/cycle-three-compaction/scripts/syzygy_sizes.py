"""Measure the cyclic covariant syzygy (python-flint).

Reads the covariant formulas of data/cycle-three/CubicCovariants.lean and reports
term counts of H, Theta, J and of J^2 = 4Theta^3 + 4 s4 Theta H^4 - s6 H^6 + U*Q,
for the generic six-coefficient family and for the physical specialization.
Also checks the single weight grade, the digit bounds and the coefficient mass
used by FullSyzygy.lean.

usage: python3 syzygy_sizes.py <repo-root>      (pip install python-flint)
"""
import re, sys
import flint

root = sys.argv[1] if len(sys.argv) > 1 else '.'
src = open(root+'/data/cycle-three/CubicCovariants.lean').read()
defs = {m.group(1): re.sub(r'\^\((\d+):ℕ\)', r'**\1', re.sub(r'\((\d+):K\)', r'\1', m.group(2)))
        for m in re.finditer(r'def (\w+Coeff\d+) \(a b c d e f : K\) : K :=\s*(.*)', src)}


def covariants(a, b, c, d, e, f, x, y, z):
    env = dict(a=a, b=b, c=c, d=d, e=e, f=f)
    out = {}
    for name in ('H', 'Theta', 'J'):
        tot = 0*a
        for key, expr in defs.items():
            if key.startswith(name+'Coeff'):
                i, j, k = map(int, key[len(name)+5:])
                tot += eval(expr, {}, env)*x**i*y**j*z**k
        out[name] = tot
    return out['H'], out['Theta'], out['J']


def syzygy(label, a, b, c, d, e, f, x, y, z):
    H, Th, J = covariants(a, b, c, d, e, f, x, y, z)
    U = a*x**3+b*y**3+c*z**3+d*x**2*y+e*y**2*z+f*z**2*x
    s = a*e**2*f+b*d*f**2+c*d**2*e
    t = 27*a**2*b**2*c**2+4*a**2*c*e**3+4*a*b**2*f**3+6*a*b*c*d*e*f+4*b*c**2*d**3-d**2*e**2*f**2
    L = J**2
    R0 = 4*Th**3+4*s*Th*H**4-t*H**6
    Q, rem = divmod(L-R0, U)
    print(f'{label}: H {len(H)}, Theta {len(Th)}, J {len(J)} terms; J^2 {len(L)}, '
          f'4Θ³+4s4ΘH⁴−s6H⁶ {len(R0)}, U-quotient {len(Q)}; exact={rem == 0}')
    return H, Th, J, U, s, t


g = flint.fmpz_mpoly_ctx.get(('a', 'b', 'c', 'd', 'e', 'f', 'x', 'y'), 'lex').gens()
a, b, c, d, e, f, x, y = g
one = a**0
H, Th, J, U, s, t = syzygy('generic family (z=1)', a, b, c, d, e, f, x, y, one)

# Grade and digit bounds for the nine-variable encoding (z recovered by homogeneity, degree 18).
R = 4*Th**3+4*s*Th*H**4-t*H**6+U*(27*t*U*Th**2-108*s**3*Th*U**3+54*s*t*Th*U**2*H-72*s**2*Th*U*H**2
    -18*t*Th*H**3-432*s**4*U**4*H-297*s**2*t*U**3*H**2-108*t**2*U**2*H**3+54*s*t*U*H**4-16*s**2*H**5)
assert J**2 == R
grades = set()
for m in (J**2).to_dict():
    ea, eb, ec, ed, ee, ef, ex, ey = m; ez = 18-ex-ey
    grades.add((3*ea+2*ed+ef+ey+ez, 3*eb+ed+2*ee+ex+ez, 3*ec+ee+2*ef+ex+ey, ex+ey+ez))
print('weight grades of J^2:', grades, '(bounds: first three < 38, last < 19)')
print('max exponents of d,e,f,x,y:', [max(m[i] for m in (J**2).to_dict()) for i in range(3, 8)], '(< 19 needed)')
mass = lambda P: sum(abs(int(v)) for v in P.to_dict().values())
mJ, mT, mH, mU, ms, mt = map(mass, (J, Th, H, U, s, t))
print('expression mass bound (bits):', (mJ**2+4*mT**3+4*ms*mT*mH**4+mt*mH**6+mU*(27*mt*mU*mT**2+108*ms**3*mT*mU**3
      +54*ms*mt*mT*mU**2*mH+72*ms**2*mT*mU*mH**2+18*mt*mT*mH**3+432*ms**4*mU**4*mH+297*ms**2*mt*mU**3*mH**2
      +108*mt**2*mU**2*mH**3+54*ms*mt*mU*mH**4+16*ms**2*mH**5)).bit_length(), '(base 2^50 used)')

if '--physical' in sys.argv:  # about 40 s
    h = flint.fmpq(1, 2)
    g = flint.fmpq_mpoly_ctx.get(('r0', 'r1', 'r2', 'x', 'y'), 'lex').gens()
    r0, r1, r2, x, y = g; one = r0**0
    syzygy('physical coefficients (r_i = Y_i/X_i)', one+r0**6, one+r1**6, one+r2**6,
           h*r0*r1, h*r1*r2, h*r2*r0, x, y, one)
