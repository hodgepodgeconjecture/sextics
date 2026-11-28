"""Low-degree Lagrange formulation of the ordinary critical trace.
Checks the Hessian/source formulas; does not compute the global trace.
"""
import json
from pathlib import Path
import sympy as s
from source_root_chart import build
from line_exactness import evaluate,derivative,inverse,mv
ROOT=Path(__file__).resolve().parent
Q,J,H=build('X1');active=list(range(7))+[10]
p=13;point=[9,10,12,3,6,8,3,1,1,2,1]
jac=[[evaluate(row[j],point,p) for j in active] for row in J]
M=jac+[[1]+[0]*7];mi=inverse(M,p)
v=mv(mi,[0]*7+[1],p)
t=s.Symbol('t');z=s.symbols('z0:4');D=s.discriminant(sum(z[j]*t**j for j in range(4)),t)
vals=dict(zip(z,[point[10],point[4],point[5],point[6]]))
grad=[0]*8
for j,k in zip([7,4,5,6],range(4)):grad[j]=int(s.diff(D,z[k]).subs(vals))%p
lamdx=mv(list(map(list,zip(*mi))),grad,p);lam=lamdx[:7];dx=lamdx[7]
assert dx==4 and sum(a*b for a,b in zip(grad,v))%p==dx
hq=[[[evaluate(H[q][i][j],point,p) for j in active] for i in active] for q in range(7)]
curv=[sum(v[i]*h[i][j]*v[j] for i in range(8) for j in range(8))%p for h in hq]
acc=mv(mi,[-q for q in curv]+[0],p)
hd=[[0]*8 for _ in range(8)]
for i,k in zip([7,4,5,6],range(4)):
 for j,l in zip([7,4,5,6],range(4)):hd[i][j]=int(s.diff(D,z[k],z[l]).subs(vals))%p
hessian=sum(v[i]*hd[i][j]*v[j] for i in range(8) for j in range(8))%p
second=(hessian+sum(a*b for a,b in zip(grad,acc)))%p
lagrangian=(hessian-sum(a*b for a,b in zip(lam,curv)))%p
assert second==lagrangian
A=s.Poly((point[7]+point[2]*t)**3*(point[0]+point[1]*t)**3,t)
rhs=[int(7*A.nth(j))%p for j in range(7)]
source=mv(mi,rhs+[0],p)
dot=sum(a*b for a,b in zip(grad,source))%p
assert dot==sum(a*b for a,b in zip(lam,rhs))%p==4
# Actual degrees after fixing K,L,k1; source matrices before determinant expansion.
degrees=[max(sum(e[j] for j in active) for co,e in q if co) for q in Q]
report={'incidence_degrees':degrees,'variables':{'curve_coordinates':8,'multipliers':7,'tangent_coordinates':8,'total':23},
'equations':{'incidence':7,'stationarity':8,'tangent':7,'normalization':1,'total':23},
'maximum_equation_degree_bound':max(degrees),
'witness':{'prime':p,'point':point,'lambda':lam,'Delta_xi':dx,'Delta_xixi':second,'Delta_dot':dot},
'checks':{'source_multiplier_identity':'PASS','Hessian_second_derivative':'PASS','tangent_normalization':'PASS'},
'scope':'The witness is deliberately not critical (Delta_xi=4); it verifies the general identities. No critical-point enumeration or global trace evaluation.'}
(ROOT/'critical-lagrange.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
