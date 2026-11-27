"""Exact incidence-Jacobian and discriminant-coordinate audit."""
import json
from pathlib import Path
import sympy as s
from source_root_chart import build,Adic
from line_exactness import evaluate,inverse,mv,compile_poly,derivative,mul,power
ROOT=Path(__file__).resolve().parent

def symbolic():
    t=s.Symbol('t'); q=s.symbols('y0 y1 x0 x1 y2 z');y0,y1,x0,x1,y2,z=q
    k=s.Rational(5,7)
    H=sum(v**6 for v in q[:5])+z*z+k*(x0**3*x1*y0*y1+x0**3*z+x0*x0*y0*y1*y2*y2+x0*x1*y1**3*y2+x0*z*y1*y2+x1*z*y0*y1)
    xi,a,b,c,d,e,f,k0,k1,k2,k3=vs=s.symbols('xi a b c d e f k0 k1 k2 k3')
    sub=dict(zip(q,[xi+a*t,t,k0+b*t,k1,k2+c*t,k3+d*t+e*t*t+f*t**3]))
    P=s.Poly(s.expand(7*H.subs(sub,simultaneous=True)),t)
    eq=s.Matrix([P.nth(i) for i in range(7)])
    M=eq.jacobian(vs[:7])
    assert list(M[0,:])==[42*xi**5,0,0,0,0,0,0]
    Ps=[s.Poly(s.expand(s.diff(7*H,v).subs(sub,simultaneous=True)),t) for v in [y0,x0,y2,z]]
    W=Ps[3].as_expr()
    polys=[p.as_expr() for p in Ps[:3]]+[W,t*W,t*t*W]
    reduced=s.Matrix([[s.expand(p).coeff(t,i) for p in polys] for i in range(6)])
    assert M[1:,1:]==reduced
    # Universal column operations proving det(M6)=lc(W)^3 det(remainders).
    w0,w1,w2,w3=s.symbols('w0 w1 w2 w3',nonzero=True)
    W=w0+w1*t+w2*t*t+w3*t**3
    pp=[sum(s.Symbol(f'p{j}_{i}')*t**i for i in range(6)) for j in range(3)]
    divs=[s.div(p,W,t) for p in pp]
    quotient=s.Matrix([[s.expand(dv[0]).coeff(t,i) for dv in divs] for i in range(3)])
    original=s.Matrix([[s.expand(p).coeff(t,i) for p in pp+[W,t*W,t*t*W]] for i in range(6)])
    operation=s.eye(6);operation[3:,:3]=-quotient
    transformed=(original*operation).applyfunc(s.cancel)
    assert transformed[3:,:3]==s.zeros(3)
    remainder=s.Matrix([[s.expand(dv[1]).coeff(t,i) for dv in divs] for i in range(3)])
    assert (transformed[:3,:3]-remainder).applyfunc(s.cancel)==s.zeros(3)
    assert transformed[3:,3:].det()==w3**3
    assert operation.det()==1
    # Exact determinant cancellation in the source and the normal derivative.
    T,Nj,Ns,J,Js,Delta=s.symbols('T Nj Ns J Js Delta',nonzero=True)
    Dj=Nj/T;Ds=Ns/T
    assert s.cancel(J*Ds/(3*Dj)-J*Ns/(3*Nj))==0
    B=Js/3-J*Ds/(18*Delta)
    assert s.cancel(B/Dj-(T*Js/(3*Nj)-J*Ns/(18*Nj*Delta)))==0
    data=json.loads((ROOT/'source-deformation-exact.json').read_text())
    alpha,beta,gamma,kappa=s.symbols('alpha beta gamma k')
    loc={'alpha':alpha,'beta':beta,'gamma':gamma,'k':kappa}
    Ds7=s.sympify(data['leading_discriminant_deformation'],locals=loc)
    D4=-27*alpha**2*gamma**2
    b04=s.factor((-s.I/6)*Ds7/(3*s.diff(D4,gamma)))
    R3=s.factor(-alpha**3/54-s.diff(b04,gamma))
    C3=s.sympify(data['compact_C0_lambda3'],locals=loc)
    assert s.cancel(C3/9-(s.diff(b04,gamma)-b04/(3*gamma))-R3)==0
    assert s.diff(R3,gamma,2)==0 and R3!=0
    correction=s.Rational(3,2)*R3.subs(gamma,0)*gamma+s.Rational(3,5)*s.diff(R3,gamma)*gamma**2
    assert s.cancel(s.diff(correction,gamma)-correction/(3*gamma)-R3)==0
    return {'incidence_matrix_structure':'PASS','universal_remainder_determinant':'PASS',
            'determinant_cancellation':'PASS','leading_normal_current':'PASS',
            'leading_pole_current':str(b04),'leading_regular_residual':str(R3),
            'leading_correction':str(s.factor(correction))}

def det_mod(a,mod):
    a=[list(r) for r in a];out=1
    for j in range(len(a)):
        h=next(i for i in range(j,len(a)) if s.gcd(a[i][j],mod)==1)
        if h!=j:a[h],a[j]=a[j],a[h];out=-out
        v=a[j][j];out=out*v%mod
        for i in range(j+1,len(a)):
            c=a[i][j]*pow(v,-1,mod)%mod
            a[i]=[(x-c*y)%mod for x,y in zip(a[i],a[j])]
    return out%mod

def audit(p,n=32):
    Adic.p=p;Adic.cap=2*n;mod=p**n
    Q,D,_=build('X1');ri=next(i for i in range(p) if (i*i+1)%p==0)
    point=[ri,0,0,ri,0,0,0,p,1,2*p,3*p]
    for _ in range(n.bit_length()+2):
        jac=[[evaluate(v,point,mod) for v in row[:7]] for row in D]
        corr=mv(inverse(jac,mod),[evaluate(q,point,mod) for q in Q],mod)
        point[:7]=[(x-y)%mod for x,y in zip(point[:7],corr)]
    assert all(evaluate(q,point,mod)==0 for q in Q)
    jac=[[evaluate(v,point,mod) for v in row] for row in D]
    T=det_mod([r[:7] for r in jac],mod)
    adj=[[T*x%mod for x in row] for row in inverse([r[:7] for r in jac],mod)]
    K=mv(adj,[-r[10] for r in jac],mod)+[0,0,0,T]
    t=s.Symbol('t');a,b,c,d=s.symbols('a b c d');zvars=[a,b,c,d]
    dc=compile_poly(s.discriminant(a*t**3+b*t*t+c*t+d,t),zvars)
    Z=[point[6],point[5],point[4],point[10]]
    grad=[evaluate(derivative(dc,j),Z,mod) for j in range(4)]
    Delta=evaluate(dc,Z,mod)
    A=mul(power([point[7],point[2]],3,mod),power([point[0],point[1]],3,mod),mod)
    A+=[0]*(7-len(A));S=mv(adj,[7*v for v in A],mod)+[0]*4
    assert all(sum(x*y for x,y in zip(row,K))%mod==0 for row in jac)
    assert all((sum(x*y for x,y in zip(row,S))-7*T*v)%mod==0 for row,v in zip(jac,A))
    Nj=sum(x*y for x,y in zip(grad,[K[6],K[5],K[4],K[10]]))%mod
    Ns=sum(x*y for x,y in zip(grad,[S[6],S[5],S[4],S[10]]))%mod
    # The Jacobian of (Q_0,...,Q_6, Delta) in (dependent variables,k3)
    # has determinant Nj by the block determinant formula.
    grow=[0]*11
    for j,v in zip([6,5,4,10],grad):grow[j]=v
    # Avoid dividing by Nj: expand the final row after eliminating its
    # first seven entries with the invertible incidence block.
    lower=mv([grow[:7]],mv(inverse([r[:7] for r in jac],mod),[r[10] for r in jac],mod),mod)[0]
    assert T*(grow[10]-lower)%mod==Nj
    vals={name:Adic(v,n) for name,v in [('T',T),('Nj',Nj),('Ns',Ns),('Delta',Delta)]}
    assert [vals[q].val() for q in ['T','Nj','Ns','Delta']]==[0,3,7,4]
    xi=Adic(point[0],n);J=1/(6*xi**5)
    old=J*(vals['Ns']/vals['T'])/(3*(vals['Nj']/vals['T']))
    new=J*vals['Ns']/(3*vals['Nj'])
    assert old==new
    return {'prime':p,'input_precision':n,'checks':{'incidence':'PASS','tangent_cofactor_field':'PASS',
        'source_cofactor_field':'PASS','new_coordinate_Jacobian':'PASS','pole_current_cancellation':'PASS'},
        'valuations':{k:v.val() for k,v in vals.items()},'pole_current_valuation':new.val()}

if __name__=='__main__':
    report={'symbolic_checks':symbolic(),'samples':[audit(p) for p in [101,109]],
      'scope':'Coordinate-chart and exact cancellation results only. No global singular-locus classification, pole bound, rational primitive or target certificate.'}
    (ROOT/'discriminant-chart.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
