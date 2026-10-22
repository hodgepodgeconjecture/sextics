"""Normal forms of the cyclic covariants modulo U (chart z=1): eliminate c, which U contains linearly.

usage: python3 reduce_mod_u.py <repo-root>"""
import re, sys, flint
src = open(sys.argv[1]+'/data/cycle-three/CubicCovariants.lean').read()
defs = {m.group(1): re.sub(r'\^\((\d+):ℕ\)', r'**\1', re.sub(r'\((\d+):K\)', r'\1', m.group(2)))
        for m in re.finditer(r'def (\w+Coeff\d+) \(a b c d e f : K\) : K :=\s*(.*)', src)}
ctx = flint.fmpz_mpoly_ctx.get(('a','b','c','d','e','f','x','y'), 'lex'); a,b,c,d,e,f,x,y = ctx.gens()
def cov(n, cc):
    env = dict(a=a,b=b,c=cc,d=d,e=e,f=f); t = 0*a
    for k, s in defs.items():
        if k.startswith(n+'Coeff'):
            i,j,_ = map(int, k[len(n)+5:]); t += eval(s, {}, env)*x**i*y**j
    return t
cstar = -(a*x**3+b*y**3+d*x**2*y+e*y**2+f*x)
for name in ('H','Theta','J'):
    full, red = cov(name, c), cov(name, cstar)
    print(f'{name:6s} terms: original {len(full):6d}   with c eliminated mod U {len(red):7d}')
s4 = lambda cc: a*e**2*f+b*d*f**2+cc*d**2*e
s6 = lambda cc: 27*a**2*b**2*cc**2+4*a**2*cc*e**3+4*a*b**2*f**3+6*a*b*cc*d*e*f+4*b*cc**2*d**3-d**2*e**2*f**2
H, Th, J = cov('H', cstar), cov('Theta', cstar), cov('J', cstar)
L, R = J**2, 4*Th**3+4*s4(cstar)*Th*H**4-s6(cstar)*H**6
print('reduced syzygy holds exactly (no U-multiple):', L == R, '; J^2 terms', len(L))
