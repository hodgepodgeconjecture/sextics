"""Read msolve -P 1 output over Q and factor the eliminating polynomial (Galois orbits)."""
import sys, re, flint
txt = open(sys.argv[1]).read()
# Output: [0, [0, nvars, deg, [vars], [linform], [1, [[deg, [coeffs of elim]], ...
m = re.search(r'\[1,\s*\[\[(\d+),\s*\[([^\]]*)\]', txt)
deg = int(m.group(1)); cs = [int(c) for c in m.group(2).replace('\n','').split(',')]
P = flint.fmpz_poly(cs)
fac = P.factor()
degs = sorted((f.degree() for f, e in fac[1] for _ in range(e)))
print('eliminant degree', P.degree(), '-> factor degrees', degs)
