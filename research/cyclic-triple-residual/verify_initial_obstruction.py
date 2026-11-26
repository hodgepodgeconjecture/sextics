"""Characteristic-zero first-jet proof, independent of the p-adic evaluator.

The scalar r is the endpoint Y0 coordinate, not a covering coordinate.
All calculations use exact Q(i) arithmetic. The coupling k stays symbolic.
"""
import json
from pathlib import Path
import sympy as s

t, r, k, alpha, beta, gamma = s.symbols('t r k alpha beta gamma', nonzero=True)
I=s.I
# At k0=k2=k3=0, k1=1, H= r^6+1 = epsilon.
baseline=s.Matrix([r,t,0,1,I*t,0])
q=s.symbols('y0 y1 x0 x1 y2 z')
y0,y1,x0,x1,y2,z=q
H=sum(v**6 for v in q[:5])+z*z+k*(x0**3*x1*y0*y1+x0**3*z+
    x0*x0*y0*y1*y2*y2+x0*x1*y1**3*y2+x0*z*y1*y2+x1*z*y0*y1)
sub=dict(zip(q,baseline))
assert s.expand(H.subs(sub,simultaneous=True))==r**6+1

# First lambda jet at endpoints (k0,k1,k2,k3)=(alpha*lambda,1,beta*lambda,gamma*lambda).
X=alpha-6*beta/k*t
Z=gamma-alpha*I/r*t**3
jet=s.Matrix([-k*gamma/(6*r**4)*t,0,X,0,beta,Z])
gradient=s.Matrix([s.diff(H,v).subs(sub,simultaneous=True) for v in q])
assert s.simplify(gradient.dot(jet))==0

# Seven unknown coefficients: xi,a,b,c,d,e,f; equations [t^0,...,t^6]H.
variation=s.Matrix([[1,t,0,0,0,0,0], [0,0,0,0,0,0,0],
                   [0,0,t,0,0,0,0], [0,0,0,0,0,0,0],
                   [0,0,0,t,0,0,0], [0,0,0,0,t,t*t,t**3]])
row=(gradient.T*variation)
incidence=s.Matrix([[s.expand(row[j]).coeff(t,n) for j in range(7)] for n in range(7)])
assert s.simplify(incidence.det()) != 0

# After removing the Y1/time and X1/k1 rows/columns, the evaluation determinant
# has rows (Y0,X0,Y2,Z), columns (epsilon,k0,k2,k3), all at lambda=0.
M=s.Matrix([[1/(6*r**5),0,0,-k/(6*r**4)*t],
            [0,1,-6/k*t,0], [0,0,1,0], [0,-I/r*t**3,0,1]])
J=s.factor(M.det()); assert J==1/(6*r**5)
normal=lambda v:s.diff(v,r)/(6*r**5)
Jdot=normal(J)
A3=X**3*r**3
AJ3=s.expand(A3*J)
AJdot3=s.expand(normal(A3)*J+A3*Jdot)
# Z is lambda Z_1 + ..., A is lambda^3 A_3 + ... .
R4=s.factor((3*Z*AJdot3-2*AJ3*normal(Z)).subs(r,I))
expected=(alpha*t**3/s.Integer(9)-gamma/s.Integer(6))*X**3
assert s.simplify(R4-expected)==0
Z1=s.simplify(Z.subs(r,I))
assert Z1==gamma-alpha*t**3
obstruction=s.expand(R4).coeff(t,4)
assert s.simplify(obstruction+2*alpha**3*beta/k)==0

# A polynomial s of degree >=5 has L(s) of degree deg(s)+2 >=7.
# For degree <=4, the coefficient of t^4 in L(s) is identically zero.
L=lambda v:s.expand(3*Z1*s.diff(v,t)-2*s.diff(Z1,t)*v)
assert all(L(t**n).coeff(t,4)==0 for n in range(5))

# Numeric rational instance, plus exact left-nullspace witness.
spec={alpha:1,beta:2,gamma:3,k:s.Rational(5,7)}
RR=s.expand(R4.subs(spec)); ZZ=Z1.subs(spec)
LL=lambda v:s.expand(3*ZZ*s.diff(v,t)-2*s.diff(ZZ,t)*v)
matrix=s.Matrix([[LL(t**j).coeff(t,i) for j in range(5)] for i in range(7)])
rhs=s.Matrix([RR.coeff(t,i) for i in range(7)])
dual=s.Matrix([0,0,0,0,1,0,0])
assert dual.T*matrix==s.zeros(1,5)
assert (dual.T*rhs)[0]==-s.Rational(28,5)
assert matrix.rank()==5 and matrix.row_join(rhs).rank()==6

# Complete two-obstruction normal form, valid for every squarefree cubic.
za,zb,zc,zd=s.symbols('za zb zc zd')
ZZ=za*t**3+zb*t*t+zc*t+zd
LL=lambda v:s.expand(3*ZZ*s.diff(v,t)-2*s.diff(ZZ,t)*v)
cols=[LL(1),LL(t),LL(t*t),ZZ,s.expand(t*ZZ)]
normal_matrix=s.Matrix([[v.coeff(t,n) for v in cols] for n in range(5)])
assert s.expand(normal_matrix.det()-8*za*s.discriminant(ZZ,t))==0
# Terms of degree >=5 reduce uniquely by L(t^n), n>=3. The five columns
# above then give the unique residual coefficients of Z and tZ.
C0=-(5*alpha**4*k**3-216*beta**3*gamma)/(36*alpha*k**3)
C1=2*alpha**2*beta/k
S=(3*beta**2*gamma/k**2-(alpha**4*k**3+216*beta**3*gamma)/(108*alpha*k**3)*t
   +alpha**2*beta/(6*k)*t*t-4*alpha*beta**2/k**2*t**3+4*beta**3/k**3*t**4)
assert s.simplify(R4-L(S)-C0*Z1-C1*t*Z1)==0
assert s.simplify(C0.subs(beta,0))==-5*alpha**3/s.Integer(36)
# On the leading family Z1 does not depend on beta. Thus both obstructions
# have polynomial beta-antiderivatives. This is a leading-order parameter
# current only; higher jets and global rationality are NOT established.
H0=-5*alpha**3*beta/s.Integer(36)+3*beta**4*gamma/(2*alpha*k**3)
H1=alpha**2*beta**2/k
assert s.simplify(s.diff(H0,beta)-C0)==0
assert s.simplify(s.diff(H1,beta)-C1)==0
assert s.diff(Z1,beta)==0

# Explicit positive control: d(s/(3w^2)) = L(s)dt/(9w^5), w^3=Z.
w=s.symbols('w',nonzero=True); ss=s.Function('s')(t); zt=s.Function('Z')(t)
primitive=ss/(3*w*w)
deriv=s.diff(primitive,t)+s.diff(primitive,w)*s.diff(zt,t)/(3*w*w)
assert s.simplify(s.expand(9*w**5*deriv).subs(w**3,zt)-(3*zt*s.diff(ss,t)-2*s.diff(zt,t)*ss))==0

# On w^3=3-t^3, all three points at infinity have w/t tending to c, c^3=-1.
# R has degree 6. Corrections in w/t start at order (1/t)^3, so the residue
# of R dt/(9w^5) is exactly -[t^4]R/(9c^5).
cc=s.symbols('c',nonzero=True)
residue=-RR.coeff(t,4)/(9*cc**5)
assert residue==s.Rational(28,45)/cc**5
report={'coupling':'5/7','exact_jet_and_normal_checks':'PASS',
        'general_obstruction':str(obstruction),'specialized_obstruction':'-28/5',
        'matrix_rank':5,'augmented_rank':6,'residue_at_each_infinity':str(residue),
        'normal_form_determinant':'8*leading_coefficient(Z)*discriminant(Z)',
        'general_compact_coefficient':str(C0),'general_residue_coefficient':str(C1),
        'specialized_compact_coefficient':str(C0.subs(spec)),
        'specialized_residue_coefficient':str(C1.subs(spec)),
        'leading_parameter_current_check':'PASS; no global rational continuation claimed',
        'R4':str(RR),'generic_coupling_scope':'k*alpha*beta*gamma != 0',
        'scope':'Obstruction to time-only integration for this anchored weighted-line family. Parameter currents remain allowed. Not Lean-verified.'}
Path(__file__).with_name('exact-result.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
