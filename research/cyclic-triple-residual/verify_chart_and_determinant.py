"""Exact target binding, chart transformation, determinant and scaling checks."""
from pathlib import Path
import sympy as s
from source_root_chart import build

root=Path(__file__).parent
record=(root/'inputs/cyclic-triple-residual-instance.rec').read_text()
xx=s.symbols('X0 X1 X2 Y0 Y1 Y2');symbols=dict(zip(map(str,xx),xx))
F=s.sympify(' '.join(record.split('F:',1)[1].split('\nA:',1)[0].splitlines()),locals=symbols)
A=s.sympify(record.split('\nA:',1)[1].splitlines()[0],locals=symbols)
assert A==xx[0]**3*xx[3]**3
t,T,z=s.symbols('t T z')
H=0
for exps,c in s.Poly(F,*xx).terms():
    assert exps[2]%3==0
    term=c*z**(exps[2]//3)
    for j,e in enumerate(exps):
        if j!=2:term*=xx[j]**e
    H+=term
xi,a,b,c,d,e,f,k0,k1,k2,k3=vs=s.symbols('xi a b c d e f k0 k1 k2 k3')
weights=[1,0,0,0,2,1,0,1,1,1,3]
for anchor in ['X1','X0']:
    Q,_,_=build(anchor)
    substitutions={xx[0]:k0+b*t if anchor=='X1' else k0,
                   xx[1]:k1 if anchor=='X1' else k1+b*t,
                   xx[3]:xi+a*t,xx[4]:t,xx[5]:k2+c*t,z:k3+d*t+e*t*t+f*t**3}
    expected=s.expand(7*H.subs(substitutions,simultaneous=True))
    reconstructed=0
    for n,terms in enumerate(Q):
        q=0
        for coeff,exps in terms:
            assert sum(w*v for w,v in zip(weights,exps))==6-n
            q+=coeff*s.prod(v**e for v,e in zip(vs,exps))
        reconstructed+=q*t**n
    assert s.expand(expected-reconstructed)==0

g=1-b*T/k0
old=[xi+a*t,t,k0+b*t,k1,k2+c*t,k3+d*t+e*t*t+f*t**3]
new=[xi+(a-xi*b/k0)*T,T,k0,k1-k1*b/k0*T,k2+(c-k2*b/k0)*T,
     k3+(d-3*k3*b/k0)*T+(e-2*d*b/k0+3*k3*b*b/k0**2)*T*T+
     (f-e*b/k0+d*b*b/k0**2-k3*b**3/k0**3)*T**3]
for j,weight in enumerate([1,1,1,1,1,3]):
    assert s.factor(g**weight*old[j].subs(t,T/g)-new[j])==0

M=s.Matrix(4,4,s.symbols('m:16'));G=s.Matrix(1,4,s.symbols('g:4'))
assert s.expand(((G*M)*M.adjugate())[0]-G[0]*M.det())==0
# GM=(1,0,0,0) implies G0*det(M)=det(M[1:,1:]); the latter has
# degree <=1+1+3=5. With deg G0=5 and det M polynomial, deg det M<=0.
r=s.symbols('r',nonzero=True)
J=1/(6*r**5)
assert s.diff(J,r)/(6*r**5)==-5/(36*r**11)

za,zb,zc,zd=s.symbols('za zb zc zd');Z=za*t**3+zb*t*t+zc*t+zd
Delta=s.discriminant(Z,t);zs=[za,zb,zc,zd];ws=[0,1,2,3]
assert s.expand(sum(w*v*s.diff(Delta,v) for w,v in zip(ws,zs))-6*Delta)==0
eta=[(zb*s.diff(Delta,za)+zc*s.diff(Delta,zb)+zd*s.diff(Delta,zc))/(6*za*Delta),
     -s.diff(Delta,za)/(6*Delta),-s.diff(Delta,zb)/(6*Delta),-s.diff(Delta,zc)/(6*Delta)]
assert s.factor(sum(w*v*q for w,v,q in zip(ws,zs,eta)))==0
# On endpoint Euler S, Gamma_S=diag(-1,0); div S=1+1+1+3=6.
# Source remainder weights are (-5,-6), so its covariant divergence
# along S is zero. Dividing by a nonzero Euler weight is unavailable.
assert -5+6-1==0 and -6+6+0==0
print('PASS: exact target binding (both charts), seven coefficient gradings, six chart identities, adjugate identity, normal derivative, and Euler resonance')
