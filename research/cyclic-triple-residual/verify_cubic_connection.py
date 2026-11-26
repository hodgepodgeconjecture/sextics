"""Exact Gauss-Manin reduction for w^3 = a*t^3+b*t^2+c*t+d.

Basis: (dt/w^2, t*dt/w^2), allowing poles only at the three infinities.
Column convention: partial_v(omega_j) = sum_i Gamma_v[i,j] omega_i
                                + d_t(S_vj/w^2).
All identities are checked over Q(a,b,c,d), including the exact corrections.
"""
import itertools
import json
from pathlib import Path
import sympy as s

t,a,b,c,d=s.symbols('t a b c d')
variables=[a,b,c,d]
Z=a*t**3+b*t*t+c*t+d
Delta=s.discriminant(Z,t)
L=lambda q:s.expand(3*Z*s.diff(q,t)-2*s.diff(Z,t)*q)
columns=[L(1),L(t),L(t*t),Z,s.expand(t*Z)]
M=s.Matrix([[v.coeff(t,n) for v in columns] for n in range(5)])
assert s.expand(M.det()-8*a*Delta)==0
Mi=M.inv()
connection={}; corrections={}
for v in variables:
    G=s.zeros(2); primitives=[]
    for j,P in enumerate([s.Integer(1),t]):
        R=s.expand(-2*P*s.diff(Z,v))
        solution=Mi*s.Matrix([R.coeff(t,n) for n in range(5)])
        primitive=sum(s.factor(solution[n])*t**n for n in range(3))
        for i in range(2): G[i,j]=s.factor(solution[3+i]/3)
        assert s.factor(R-L(primitive)-3*Z*(G[0,j]+t*G[1,j]))==0
        primitives.append(primitive)
    assert s.factor(G[0,0]+s.diff(Delta,v)/(6*Delta))==0
    assert G[1,0]==0
    assert s.factor(G[1,1]+s.Rational(2,3)*s.diff(a,v)/a)==0
    connection[v]=G; corrections[v]=primitives

eta=[(b*s.diff(Delta,a)+c*s.diff(Delta,b)+d*s.diff(Delta,c))/(6*a*Delta),
     -s.diff(Delta,a)/(6*Delta),-s.diff(Delta,b)/(6*Delta),-s.diff(Delta,c)/(6*Delta)]
for v,e in zip(variables,eta): assert s.factor(connection[v][0,1]-e)==0
for u,v in itertools.combinations(variables,2):
    curvature=connection[v].diff(u)-connection[u].diff(v)+connection[u]*connection[v]-connection[v]*connection[u]
    assert all(s.factor(x)==0 for x in curvature)

# Scaling controls: scaling the cubic by lambda multiplies both basis
# differentials by lambda^(-2/3). Moving the time coordinate additionally
# changes its differential; that operation is separate from this control.
assert all(s.factor(x)==0 for x in sum((v*connection[v] for v in variables),s.zeros(2))+s.Rational(2,3)*s.eye(2))

out={'discriminant':str(Delta),'normal_form_determinant':'8*a*Delta',
     'basis':['dt/w^2','t*dt/w^2'],
     'connection':{str(v):[[str(x) for x in row] for row in connection[v].tolist()] for v in variables},
     'exact_time_corrections':{str(v):[str(s.factor(x)) for x in corrections[v]] for v in variables},
     'checks':{'eight_differential_identities':'PASS','six_flatness_identities':'PASS','scaling':'PASS'},
     'scope':'Universal rational differential identities, not a source-current solution or Lean proof.'}
Path(__file__).with_name('connection.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'checks':out['checks'],'Gamma00':'-dlog(Delta)/6','Gamma11':'-2*dlog(a)/3',
                  'Gamma10':0,'Gamma01':'[(b*Delta_a+c*Delta_b+d*Delta_c)*da/a-Delta_a*db-Delta_b*dc-Delta_c*dd]/(6*Delta)'},indent=2))
