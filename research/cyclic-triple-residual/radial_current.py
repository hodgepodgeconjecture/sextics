"""Construct a formal radial current; no rationality or certificate claim."""
import json
from pathlib import Path
import sympy as s
from source_root_chart import build
from line_exactness import compile_poly, inverse
ROOT=Path(__file__).resolve().parent

def pole_removal_identity():
    J,Js,Jj,D,Ds,Dj,Dsj,Djj=s.symbols('J Js Jj D Ds Dj Dsj Djj')
    b=J*Ds/(3*Dj)
    db=s.diff(b,J)*Jj+s.diff(b,Ds)*Dsj+s.diff(b,Dj)*Djj
    source=Js/3-J*Ds/(18*D)
    assert s.cancel(source-(db-b*Dj/(6*D))-(Js/3-db))==0
    C,Cj=s.symbols('C Cj')
    assert s.cancel(D*Cj+Dj*C-D*C*Dj/(6*D)-(D*Cj+5*Dj*C/6))==0

def audit(p,N,direction):
    def const(x): return [int(x)%p]+[0]*(N-1)
    def add(a,b): return [(x+y)%p for x,y in zip(a,b)]
    def scale(a,k): return [x*k%p for x in a]
    def mul(a,b):
        c=[0]*N
        for i,x in enumerate(a):
            if x:
                for j in range(N-i): c[i+j]=(c[i+j]+x*b[j])%p
        return c
    def inv(a):
        assert a[0]
        b=const(pow(a[0],-1,p))
        for n in range(1,N): b[n]=-b[0]*sum(a[j]*b[n-j] for j in range(1,n+1))%p
        return b
    def power(a,n):
        out=const(1)
        for _ in range(n): out=mul(out,a)
        return out
    def ev_many(polys,point):
        powers=[[const(1),v] for v in point]
        for row in powers:
            for _ in range(5): row.append(mul(row[-1],row[1]))
        out=[]
        for terms in polys:
            acc=const(0)
            for co,exps in terms:
                term=const(co)
                for j,e in enumerate(exps):
                    if e: term=mul(term,powers[j][e])
                acc=add(acc,term)
            out.append(acc)
        return out
    Q,D,_=build('X1')
    ri=next(i for i in range(p) if (i*i+1)%p==0)
    point=list(map(const,[ri,0,0,ri,0,0,0,0,1,0,0]))
    for j,a in zip([7,9,10],direction): point[j][1]=a%p
    flat=[v for row in D for v in row[:7]]
    J0=ev_many(flat,point)
    ji=inverse([[J0[7*i+j][0] for j in range(7)] for i in range(7)],p)
    # Constant-Jacobian Hensel iteration increases precision at least one
    # order each step. Endpoint parameters remain exactly linear.
    for _ in range(N):
        residual=ev_many(Q,point)
        for i in range(7):
            corr=[sum(ji[i][j]*residual[j][n] for j in range(7))%p for n in range(N)]
            point[i]=add(point[i],scale(corr,-1))
    assert all(not any(q) for q in ev_many(Q,point))
    jac=ev_many(flat,point)
    t=s.Symbol('t'); vs=s.symbols('xi a b c d e f k0 k1 k2 k3')
    Ap=s.Poly((vs[7]+vs[2]*t)**3*(vs[0]+vs[1]*t)**3,t)
    Avec=ev_many([compile_poly(Ap.nth(j),vs) for j in range(7)],point)
    U=[const(0) for _ in range(7)]
    for n in range(N):
        rhs=[(7*Avec[i][n]-sum(jac[7*i+j][h]*U[j][n-h]
                for j in range(7) for h in range(1,n+1)))%p for i in range(7)]
        for i in range(7): U[i][n]=sum(ji[i][j]*rhs[j] for j in range(7))%p
    for i in range(7):
        got=const(0)
        for j in range(7): got=add(got,mul(jac[7*i+j],U[j]))
        assert got==scale(Avec[i],7)
    a,b,c,d=s.symbols('a b c d'); variables=[a,b,c,d]
    disc=s.discriminant(a*t**3+b*t*t+c*t+d,t)
    dd=ev_many([compile_poly(disc,variables)]+
               [compile_poly(s.diff(disc,v),variables) for v in variables],
               [point[6],point[5],point[4],point[10]])
    Delta=dd[0]; Ds=const(0)
    for gradient,variation in zip(dd[1:],[U[6],U[5],U[4],const(0)]):
        Ds=add(Ds,mul(gradient,variation))
    assert Delta[:4]==[0]*4 and Delta[4]
    assert Ds[:7]==[0]*7
    unit=Delta[4:]+[0]*4
    ratio=mul(Ds[4:]+[0]*4,inv(unit))
    J=scale(inv(power(point[0],5)),pow(6,-1,p))
    Js=scale(mul(power(point[7],3),inv(power(point[0],8))),-pow(18,-1,p))
    B=add(scale(Js,pow(3,-1,p)),scale(mul(J,ratio),-pow(18,-1,p)))
    # Only N-4 terms are known after dividing by the order-four Delta.
    precision=N-4
    logunit=mul([n*unit[n]%p for n in range(N)],inv(unit))
    coeff=scale(logunit,-pow(6,-1,p)); coeff[0]=7*pow(3,-1,p)%p
    f=const(0)
    for n in range(precision):
        denominator=(n+coeff[0])%p
        assert denominator
        f[n]=(B[n]-sum(coeff[j]*f[n-j] for j in range(1,n+1)))*pow(denominator,-1,p)%p
    residual=add([n*f[n]%p for n in range(N)],mul(coeff,f))
    assert residual[:precision]==B[:precision]
    assert f[:3]==[0]*3 and f[3]
    exact=json.loads((ROOT/'source-deformation-exact.json').read_text())
    alpha,beta,gamma,k=s.symbols('alpha beta gamma k')
    lead=s.sympify(exact['compact_C0_lambda3'],locals={'alpha':alpha,'beta':beta,'gamma':gamma,'k':k})/48
    val=s.cancel(lead.subs({alpha:direction[0],beta:direction[1],gamma:direction[2],k:s.Rational(5,7),s.I:ri}))
    num,den=s.fraction(val)
    assert f[3]==int(num)*pow(int(den),-1,p)%p
    return {'prime':p,'direction':direction,'incidence_precision':N,
            'radial_equation_precision':precision,'source_coefficients':B[:precision],
            'radial_f_coefficients':f[:precision],
            'checks':{'incidence':'PASS','source_variation':'PASS','radial_current_equation':'PASS',
                      'independent_characteristic_zero_leading_current':'PASS'}}

def leading_identity():
    data=json.loads((ROOT/'source-deformation-exact.json').read_text())
    x,y,z,k=s.symbols('alpha beta gamma k')
    C=s.sympify(data['compact_C0_lambda3'],locals={'alpha':x,'beta':y,'gamma':z,'k':k})
    f=C/48; Delta=-27*x*x*z*z
    divergence=sum(s.diff(q*f,q)-q*f*s.diff(Delta,q)/(6*Delta) for q in [x,y,z])
    assert s.cancel(divergence-C/9)==0
    return str(s.factor(f))

if __name__=='__main__':
    pole_removal_identity()
    leading=leading_identity()
    runs=[audit(p,24,d) for p in [101,109] for d in [(1,2,3),(2,1,4)]]
    report={'rational_discriminant_pole_removal':'PASS',
            'leading_radial_coefficient_f3':leading,'runs':runs,
            'scope':'Formal radial current only. Finite-field jets do not prove rationality, descent, convergence, or HasCertificate.'}
    (ROOT/'radial-current.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'leading_identity':'PASS','runs':[{'prime':r['prime'],'direction':r['direction'],
          'precision':r['radial_equation_precision'],'checks':r['checks']} for r in runs]},indent=2))
