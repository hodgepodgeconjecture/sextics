"""Frobenius cycle types of the 240 minimal sections (factor the mod-p eliminant)."""
import re, subprocess, collections, sys, flint
c1, c2 = sys.argv[1:3]
primes = [p for p in range(1009, 40000) if flint.fmpz(p).is_prime()][:60]
types = collections.Counter()
for p in primes:
    subprocess.run(['python3', 'gen_sections.py', c1, c2, str(p), 's.ms'], check=True, capture_output=True)
    subprocess.run(['msolve', '-f', 's.ms', '-o', 's.out', '-P', '1'], check=True, capture_output=True)
    txt = open('s.out').read()
    m = re.search(r'\[1,\s*\[\[(\d+),\s*\[([^\]]*)\]', txt)
    if not m or int(m.group(1)) != 240: continue
    P = flint.nmod_poly([int(c) for c in m.group(2).replace('\n', '').split(',')], p)
    fac = P.factor()[1]
    if any(e > 1 for _, e in fac): continue          # skip ramified primes
    types[tuple(sorted((f.degree() for f, _ in fac), reverse=True))] += 1
for t, n in types.most_common():
    c = collections.Counter(t)
    print(f'{n:3d}x  ' + ' '.join(f'{d}^{k}' for d, k in sorted(c.items(), reverse=True)))
