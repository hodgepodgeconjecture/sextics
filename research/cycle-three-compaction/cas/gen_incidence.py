"""Write the 13x13 cycle-three Jacobian incidence system in msolve format."""
import sys, flint
c1, c2, R0, X0 = map(int, sys.argv[1:5]); char = int(sys.argv[5]); out = sys.argv[6]
names = [f'v{i}' for i in range(13)] + ['T']
ctx = flint.fmpz_mpoly_ctx.get(tuple(names), 'degrevlex'); g = ctx.gens(); v, T = g[:13], g[13]
one = T**0
def a8(c1, c2, al, r, t): return r**2*t*(c1**2*c2**3*(al*r+t)+c1*c2**2*(r+c1**6*t)+c1**3*c2*(r+c2**6*t))
def b256(c1, c2, al, r, t):
    return (-1728*(al*r+t)**2*(r+c1**6*t)**2*(r+c2**6*t)**2
            -32*r**2*t*(c1**3*c2**3*(al*r+t)**2*(r+c2**6*t)+c2**3*(al*r+t)*(r+c1**6*t)**2+c1**3*(r+c1**6*t)*(r+c2**6*t)**2)
            -48*r**2*t*c1**2*c2**2*(al*r+t)*(r+c1**6*t)*(r+c2**6*t)+r**4*t**2*c1**4*c2**4)
RR = R0*one + v[0]*T + v[1]*T**2
XX = X0*one + v[2]*T + v[3]*T**2 + v[4]*T**3 + v[5]*T**4
YY = sum(v[6+k]*T**k for k in range(7))
E = 256*YY**2 - 256*XX**3 - 32*a8(c1, c2, 1, RR, T)*XX - b256(c1, c2, 1, RR, T)
coeffs = {}
for exps, cf in E.to_dict().items():
    k = exps[13]; mono = exps[:13]
    coeffs.setdefault(k, {})[mono+(0,)] = cf
polys = []
for k in sorted(coeffs):
    P = ctx.from_dict(coeffs[k])
    polys.append(str(P).replace('^', '^'))
with open(out, 'w') as f:
    f.write(','.join(names[:13]) + '\n' + str(char) + '\n' + ',\n'.join(polys) + '\n')
print(len(polys), 'equations; degrees', [ctx.from_dict(coeffs[k]).total_degree() for k in sorted(coeffs)])
