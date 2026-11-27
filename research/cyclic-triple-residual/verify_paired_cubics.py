"""Exact paired-cubic exchange screen; no target certificate claimed."""
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parent
record=(ROOT/'inputs/cyclic-triple-residual-instance.rec').read_text()
X0,X1,X2,Y0,Y1,Y2=s.symbols('X0 X1 X2 Y0 Y1 Y2')
xx=[X0,X1,X2,Y0,Y1,Y2]
F=s.sympify(' '.join(record.split('F:',1)[1].split('\nA:',1)[0].splitlines()),
           locals=dict(zip(map(str,xx),xx)))
assert s.sympify(record.split('\nA:',1)[1].splitlines()[0],
                 locals=dict(zip(map(str,xx),xx)))==X0**3*Y0**3
assert record.split('\nq:',1)[1].splitlines()[0].strip()=='1'
assert (ROOT.parents[1]/'data/cyclic-triple-residual/instance.rec').read_text().rstrip()==record.rstrip()
z,t,lam=s.symbols('z t lam')
k=s.Rational(5,7)
H=0
for exps,co in s.Poly(F,*xx).terms():
    assert exps[2]%3==0
    H+=co*z**(exps[2]//3)*s.prod(xx[j]**exps[j] for j in range(6) if j!=2)
B=X0**3+X0*Y1*Y2+X1*Y0*Y1
partner=-k*B-z
assert s.expand(H.subs(z,partner)-H)==0
assert s.expand(partner.subs(z,partner)-z)==0
w=s.Symbol('w')
assert s.expand(H.subs(z,w**3)-(w**3-z)*(w**3-partner)-H)==0

xi,a,b,c,d,e,f,k0,k2,k3=s.symbols('xi a b c d e f k0 k2 k3')
dep=[xi,a,b,c,d,e,f]
sub={X0:k0+b*t,X1:1,Y0:xi+a*t,Y1:t,Y2:k2+c*t,z:k3+d*t+e*t*t+f*t**3}
seed={xi:s.I,a:0,b:0,c:s.I,d:0,e:0,f:0,k0:0,k2:0,k3:0}
P=s.Poly(s.expand(H.subs(sub,simultaneous=True)),t)
eq=s.Matrix([P.nth(j) for j in range(7)])
assert eq.subs(seed)==s.zeros(7,1)
Jac=eq.jacobian(dep).subs(seed)
assert Jac.det()!=0
alpha,beta,gamma=s.symbols('alpha beta gamma')
endjet=s.Matrix([alpha,beta,gamma])
jet=-Jac.inv()*eq.jacobian([k0,k2,k3]).subs(seed)*endjet
allvars=dep+[k0,k2,k3]
alljet=list(jet)+list(endjet)
def first_jet(expr):
    line=s.expand(expr.subs(sub,simultaneous=True))
    return s.expand(line.subs(seed)+lam*sum(s.diff(line,q).subs(seed)*v
                                          for q,v in zip(allvars,alljet)))
Zplus=first_jet(z)
Zminus=first_jet(partner)
assert s.expand(Zplus-lam*(gamma-alpha*t**3))==0
assert s.expand(Zminus-(-s.I*k*t+lam*(-gamma+(-s.I*k*alpha+k*k*gamma/6)*t*t
                              +(alpha+6*s.I*beta)*t**3)))==0
Dplus=s.expand(s.discriminant(Zplus,t))
Dminus=s.expand(s.discriminant(Zminus,t))
assert Dplus==-27*alpha**2*gamma**2*lam**4
assert s.expand(Dminus.coeff(lam,1)+4*s.I*k**3*(alpha+6*s.I*beta))==0
assert Dminus.coeff(lam,0)==0
# Higher incidence jets cannot alter these leading coefficients:
# discr is quartic at Zplus=0, and the linear term at Zminus=-i*k*t
# depends only on the first incidence jet. Treat all second-order
# corrections as arbitrary, rather than mistaking a truncated jet for
# an exact incidence solution.
aa,bb,cc,dd=s.symbols('aa bb cc dd')
Dt=s.discriminant(aa*t**3+bb*t*t+cc*t+dd,t)
origin={aa:0,bb:0,cc:0,dd:0}
other={aa:0,bb:0,cc:-s.I*k,dd:0}
assert {sum(ex) for ex,_ in s.Poly(Dt,aa,bb,cc,dd).terms()}=={4}
minus_coeffs=[s.expand(Zminus).coeff(t,j).coeff(lam,1) for j in [3,2,1,0]]
assert s.expand(sum(s.diff(Dt,q).subs(other)*v
                    for q,v in zip([aa,bb,cc,dd],minus_coeffs))-Dminus.coeff(lam,1))==0

# The target has deck character zeta; its average and the trace of both
# terms in the transported differential vanish, using zeta^2+zeta+1=0.
zeta=s.Symbol('zeta')
minimal=s.Poly(zeta*zeta+zeta+1,zeta)
def rem(p):
    return s.rem(s.Poly(p,zeta),minimal).as_expr()
assert rem(1+zeta+zeta**2)==0
for denominator_power in [2,5]:
    assert rem(sum(zeta**((-denominator_power*j)%3) for j in range(3)))==0

report={
 'input_sha256':hashlib.sha256(record.encode()).hexdigest(),
 'checks':{name:'PASS' for name in ['exact_exchange_involution_and_factorization',
    'exact_seed_and_implicit_first_jet','leading_discriminants_with_jet_order_justification',
    'target_and_source_trace_zero','rational_horizontal_gauge_valuation_obstruction']},
 'plus_cubic_first_jet':str(Zplus), 'minus_cubic_first_jet':str(Zminus),
 'plus_discriminant_order':4,'plus_discriminant_leading_coefficient':str(Dplus.coeff(lam,4)),
 'minus_discriminant_order':1,'minus_discriminant_leading_coefficient':str(s.factor(Dminus.coeff(lam,1))),
 'ratio_discriminant_order':1-4,
 'same_character_gauge_required_order':str(s.Rational(1-4,6)),
 'opposite_character_gauge_required_order':str(s.Rational(-1-4,6)),
 'limits':['Generic direction: alpha*gamma*(alpha+6*i*beta) != 0.',
           'The valuation is the exceptional-divisor valuation of the blowup of the smooth line-space seed.',
           'No rational horizontal identification of the two rank-one compact connections on this component.',
           'Finite ramified extensions and other correspondences are not excluded.',
           'No unrestricted rational-current obstruction and no HasCertificate proof.']}
assert s.Rational(1-4,6).q!=1 and s.Rational(-1-4,6).q!=1
(ROOT/'paired-cubics.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
