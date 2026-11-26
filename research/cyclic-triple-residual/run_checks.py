"""Replay the checkpoint and compare its independent exact and adic routes."""
import json
from pathlib import Path
import subprocess
import sys
import sympy as sp

ROOT=Path(__file__).resolve().parent
commands=[['verify_initial_obstruction.py'],['verify_cubic_connection.py'],
          ['verify_chart_and_determinant.py'],
          ['line_exactness.py','--prime','101','--precision','8','--output','result.json'],
          ['line_exactness.py','--prime','109','--precision','10','--output','result-109.json'],
          ['source_root_chart.py','--prime','101','--precision','60','--output','source-root-101.json'],
          ['source_root_chart.py','--prime','109','--precision','40','--output','source-root-109.json']]
for command in commands:
    run=subprocess.run([sys.executable,*command],cwd=ROOT,text=True,capture_output=True)
    if run.returncode:
        print(run.stdout);print(run.stderr,file=sys.stderr);raise SystemExit(run.returncode)
    print('PASS: '+' '.join(command),flush=True)

t=sp.symbols('t')
R=sp.expand((t**3/sp.Integer(9)-sp.Rational(1,2))*(1-sp.Rational(84,5)*t)**3)
for name in ['result.json','result-109.json']:
    report=json.loads((ROOT/name).read_text());p=report['prime']
    def mod(q):
        n,d=sp.fraction(q);return int(n)*pow(int(d),-1,p)%p
    assert all(v==mod(R.coeff(t,j)) for j,v in enumerate(report['normalized_R_mod_p']))
    assert report['normalized_remainders_mod_p']==[mod(sp.Rational(1777487,4500)),mod(sp.Rational(28,5))]
    assert (report['rank'],report['augmented_rank'])==(8,9)
    assert report['full_normal_form_reconstruction']
for name in ['source-root-101.json','source-root-109.json']:
    report=json.loads((ROOT/name).read_text())
    for component in ['compact_remainder','residue_remainder']:
        assert report[component]['valuation']==2
        assert report[component]['absolute_precision']>2
        assert not report[component]['zero_to_precision']
print('PASS: independent characteristic-zero / adic agreement, and source-root nonzero remainders')
print('All research checks passed. No Lean verification or target certificate claimed.')
