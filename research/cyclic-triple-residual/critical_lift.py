"""Exact p-adic lift of an ordinary discriminant critical point.
A local summand, not the global trace or an exactness verdict.
"""
import json
from pathlib import Path
import sympy as s
from algebraic_curve_model import model
from source_root_chart import build
from line_exactness import compile_poly,derivative,evaluate,inverse,mv
ROOT=Path(__file__).resolve().parent
model() # Bind the quotient construction to the protected target.
raw,_,_=build('X1')
x=s.symbols('xi a b c d e f z');lam=s.symbols('l0:7');vs=x+lam
xi,a,b,c,d,e,f,z=x;t=s.Symbol('t')
old=[xi,a,b,c,d,e,f,s.Integer(1),s.Integer(1),s.Integer(2),z]
Q=[s.expand(sum(co*s.prod(v**n for v,n in zip(old,ex)) for co,ex in q)) for q in raw]
Delta=s.discriminant(z+d*t+e*t*t+f*t**3,t)
J=s.Matrix(Q).jacobian(x)
station=[s.diff(Delta,v)-sum(lam[j]*J[j,i] for j in range(7)) for i,v in enumerate(x)]
polys=[compile_poly(q,vs) for q in Q+station]
jac=[[derivative(q,j) for j in range(15)] for q in polys]
p=13;n=32;mod=p**n
seed=[4,6,2,12,6,9,3,1,3,0,10,1,11,9,3];point=seed[:]
assert all(evaluate(q,point,p)==0 for q in polys)
seedjac=[[evaluate(q,point,p) for q in row] for row in jac]
seed_det=int(s.Matrix(seedjac).det())%p;assert seed_det
for _ in range(n.bit_length()+1):
 mat=[[evaluate(q,point,mod) for q in row] for row in jac]
 correction=mv(inverse(mat,mod),[evaluate(q,point,mod) for q in polys],mod)
 point=[(u-v)%mod for u,v in zip(point,correction)]
assert all(evaluate(q,point,mod)==0 for q in polys)
def ev(expr):return evaluate(compile_poly(s.expand(expr),vs),point,mod)
jnum=[[ev(q) for q in row] for row in J.tolist()]
M=jnum+[[1]+[0]*7];v=mv(inverse(M,mod),[0]*7+[1],mod)
H=s.hessian(Delta,x)-sum((lam[j]*s.hessian(Q[j],x) for j in range(7)),s.zeros(8))
B=sum(v[i]*ev(H[i,j])*v[j] for i in range(8) for j in range(8))%mod
A=s.Poly((1+b*t)**3*(xi+a*t)**3,t)
dot=sum(point[8+j]*ev(7*A.nth(j)) for j in range(7))%mod
# Use 7*w to keep polynomial evaluation integral, then invert 7.
w=ev(14*z+5)*pow(7,-1,mod)%mod;delta=ev(Delta)
term=2*dot*pow(w*w*delta*B,-1,mod)%mod
assert (B%p,dot%p)==(12,5)
report={'prime':p,'precision_exponent':n,'fixed_parameters':{'K':1,'L':2},
'variables':list(map(str,vs)),'seed':seed,'lift_mod_prime_power':point,
'critical_Jacobian_determinant_mod_p':seed_det,
'local_data_mod_p':{'Delta':delta%p,'w':w%p,'Delta_xixi':B%p,'Delta_dot':dot%p,'trace_summand':term%p},
'local_trace_summand_mod_prime_power':term,
'checks':{'critical_equations':'PASS','unit_critical_Jacobian':'PASS','Newton_lift_all_15_equations':'PASS','ordinary_chart_denominators':'PASS'},
'scope':'One exact lifted critical point at K=1,L=2. Nonzero local contribution does not imply nonzero global trace; conjugate and other critical points have not been summed.'}
(ROOT/'critical-lift.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
