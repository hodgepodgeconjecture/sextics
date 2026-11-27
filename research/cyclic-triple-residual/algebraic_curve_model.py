"""Eliminate cubic-lift coefficients and certify a generic degree bound.

The 120-sheet bound is over the endpoint curve, NOT over Delta.
No primitive or target certificate is claimed.
"""
import json
from pathlib import Path
import sympy as s
from source_root_chart import build
from line_exactness import evaluate,inverse,mv
from discriminant_chart import det_mod
ROOT=Path(__file__).resolve().parent

def monomials(n):
    return [(i,j,n-i-j) for i in range(n+1) for j in range(n-i+1)]

def source_transport():
    x,u,w,D,r,Dx,Ds=s.symbols('x u w D r Dx Ds',nonzero=True)
    wx=-12*x**5/w
    rx=r*Dx/(6*D)
    # u=k0; Ds denotes deformation at fixed xi, not fixed k3.
    old_Ds=Ds+u**3*Dx/(6*x**2)
    old_density=(u**3/(9*x**3*w)+old_Ds/(18*w*D))/r
    h=u**3/(18*x**2*w*r)
    dh=s.diff(h,x)+s.diff(h,w)*wx+s.diff(h,r)*rx
    new_density=2*u**3*x**3/(3*w**3*r)+Ds/(18*w*D*r)
    assert s.cancel(old_density+dh-new_density)==0
    # Also check the variable change directly from the original B.
    J=1/(6*x**5);Js=-u**3/(18*x**8)
    B=Js/3-J*old_Ds/(18*D)
    assert s.cancel(B*(-6*x**5/w)/r-old_density)==0
    return {'source_coordinate_change':'PASS','exact_endpoint_pole_subtraction':'PASS',
            'primitive_subtraction':'k0^3/(18*xi^2*w0*r)',
            'new_density':'2*k0^3*xi^3/(3*w0^3*r) + Delta_dot/(18*w0*Delta*r)',
            'Delta_dot_convention':'source deformation at fixed xi,k0,k2'}

def rank_witness(rows,p):
    a=[r[:] for r in rows];ids=list(range(len(a)));chosen=[];r=0
    for j in range(len(a[0])):
        q=next((i for i in range(r,len(a)) if a[i][j]%p),None)
        if q is None:continue
        a[r],a[q]=a[q],a[r];ids[r],ids[q]=ids[q],ids[r]
        inv=pow(a[r][j],-1,p);a[r]=[v*inv%p for v in a[r]]
        for i in range(r+1,len(a)):
            v=a[i][j]%p
            if v:a[i]=[(x-v*y)%p for x,y in zip(a[i],a[r])]
        chosen.append(ids[r]);r+=1
        if r==len(a):break
    return r,chosen

def model():
    xx=s.symbols('X0 X1 X2 Y0 Y1 Y2');X0,X1,X2,Y0,Y1,Y2=xx
    rec=(ROOT/'inputs/cyclic-triple-residual-instance.rec').read_text()
    assert rec.rstrip()==(ROOT.parents[1]/'data/cyclic-triple-residual/instance.rec').read_text().rstrip()
    F=s.sympify(' '.join(rec.split('F:',1)[1].split('\nA:',1)[0].splitlines()),locals=dict(zip(map(str,xx),xx)))
    A=s.sympify(rec.split('\nA:',1)[1].splitlines()[0],locals=dict(zip(map(str,xx),xx)))
    assert A==X0**3*Y0**3 and rec.split('\nq:',1)[1].splitlines()[0].strip()=='1'
    k=s.Rational(5,7);z=s.Symbol('z');H=0
    for ex,co in s.Poly(F,*xx).terms():
        assert ex[2]%3==0
        H+=co*z**(ex[2]//3)*s.prod(xx[j]**ex[j] for j in range(6) if j!=2)
    B=s.expand(H).coeff(z,1)/k;C=s.expand(H).coeff(z,0)
    assert s.expand(H-(z*z+k*B*z+C))==0
    assert s.expand(4*H-((2*z+k*B)**2-(k*k*B*B-4*C)))==0
    source_parameter=s.Symbol('source_parameter')
    assert s.expand(4*(H-source_parameter*A)-((2*z+k*B)**2-
                     (k*k*B*B-4*C+4*source_parameter*A)))==0
    t,xi,k0,k2,k3,a,b,c,h=s.symbols('t xi k0 k2 k3 a b c h')
    sub={X0:k0+b*t,X1:1,Y0:xi+a*t,Y1:h*t,Y2:k2+c*t}
    Bt=s.Poly(s.expand(B.subs(sub,simultaneous=True)),t)
    At=s.Poly(s.expand(A.subs(sub,simultaneous=True)),t)
    D=s.Poly(s.expand((k*k*B*B-4*C).subs(sub,simultaneous=True)),t)
    ds=[D.nth(i) for i in range(7)]
    w0=2*k3+k*k0**3
    endpoint=xi**6+k0**6+1+k2**6+k3*k3+k*k0**3*k3
    assert s.expand(w0*w0-ds[0]-4*endpoint)==0
    for i,di in enumerate(ds):
        assert all(sum(ex)==i for ex,co in s.Poly(di,a,b,c,h).terms() if co)
    # Universal square-root recurrence and the three residual equations.
    w=s.Symbol('w',nonzero=True);q=s.symbols('D1:7')
    d1,d2,d3,d4,d5,d6=q
    n2=4*w*w*d2-d1*d1
    n3=8*w**4*d3-4*w*w*d1*d2+d1**3
    W=w+d1*t/(2*w)+n2*t*t/(8*w**3)+n3*t**3/(16*w**5)
    delta=w*w+sum(v*t**(j+1) for j,v in enumerate(q))
    diff=s.Poly(s.expand(delta-W*W),t)
    assert all(s.cancel(diff.nth(i))==0 for i in range(4))
    residual=[64*w**6*d4-4*d1*n3-n2*n2,
              64*w**8*d5-n2*n3,256*w**10*d6-n3*n3]
    for i,den in enumerate([64*w**6,64*w**8,256*w**10]):
        assert s.cancel(diff.nth(i+4)*den-residual[i])==0
    return {'vars':(a,b,c,h),'params':(xi,k0,k2,k3),'D':ds,'B':[Bt.nth(i) for i in range(4)],
            'A':[At.nth(i) for i in range(7)],
            'w0':w0,'endpoint':endpoint,'residual_templates':residual,'template_vars':(w,)+q}

def local_binding(M,p,n=48):
    mod=p**n;Q,D,_=build('X1');ri=next(i for i in range(p) if (i*i+1)%p==0)
    point=[ri,0,0,ri,0,0,0,p,1,2*p,3*p]
    for _ in range(n.bit_length()+2):
        J=[[evaluate(v,point,mod) for v in row[:7]] for row in D]
        corr=mv(inverse(J,mod),[evaluate(q,point,mod) for q in Q],mod)
        point[:7]=[(x-y)%mod for x,y in zip(point[:7],corr)]
    assert all(evaluate(q,point,mod)==0 for q in Q)
    J=[[evaluate(v,point,mod) for v in row[:7]] for row in D];T=det_mod(J,mod)
    xi,k0,k2,k3=M['params'];a,b,c,h=M['vars']
    variables=[xi,k0,k2,k3,a,b,c,h]
    vals=[point[0],point[7],point[9],point[10],point[1],point[2],point[3],1]
    def ev(expr,variables=variables,vals=vals):
        out=0
        for ex,co in s.Poly(expr,*variables).terms():
            term=int(co.p)*pow(int(co.q),-1,mod)%mod
            for x,e in zip(vals,ex):term=term*pow(x,e,mod)%mod
            out=(out+term)%mod
        return out
    w=ev(M['w0']);ds=[ev(q) for q in M['D']];bs=[ev(q) for q in M['B']]
    kap=5*pow(7,-1,mod)%mod
    Ws=[(2*z+kap*b)%mod for z,b in zip([point[10],point[4],point[5],point[6]],bs)]
    n2=(4*w*w*ds[2]-ds[1]**2)%mod
    n3=(8*w**4*ds[3]-4*w*w*ds[1]*ds[2]+ds[1]**3)%mod
    assert (2*w*Ws[1]-ds[1])%mod==0
    assert (8*w**3*Ws[2]-n2)%mod==0
    assert (16*w**5*Ws[3]-n3)%mod==0
    tv=M['template_vars'];tvals=[w]+ds[1:]
    assert all(ev(f,tv,tvals)==0 for f in M['residual_templates'])
    S=[]
    for f in M['residual_templates']:
        fd=[ev(s.diff(f,q),tv,tvals) for q in tv[1:]]
        S.append([sum(u*ev(s.diff(di,x)) for u,di in zip(fd,M['D'][1:]))%mod for x in [a,b,c]])
    detS=int(s.Matrix(S).det())%mod
    assert (2**26*pow(w,21,mod)*T-42*7**6*pow(point[0],5,mod)*detS)%mod==0
    def valuation(x):
        v=0
        while x and x%p==0:x//=p;v+=1
        return v if x else n
    assert valuation(w)==1 and valuation(T)==0 and valuation(detS)==21
    return {'prime':p,'precision':n,'reconstructed_lift':'PASS','reduced_incidence':'PASS',
            'cleared_Jacobian_identity':'PASS','valuations':{'w0':1,'T':0,'det_S':21}}

def finite_test(M,p):
    xi,k0,k2,k3=M['params'];a,b,c,h=M['vars'];k=5*pow(7,-1,p)%p
    # A nonsingular point of the endpoint curve admits a characteristic-zero
    # algebraic lift with this residue, since the k3 derivative is w0 != 0.
    point=next((r,z) for r in range(1,p) for z in range(p)
       if (r**6+1+1+2**6+z*z+k*z)%p==0 and (2*z+k)%p)
    sub={xi:point[0],k0:1,k2:2,k3:point[1]}
    def reduce(poly):
        rational=s.Poly(poly,a,b,c,h)
        return s.Poly.from_dict({ex:int(co.p)*pow(int(co.q),-1,p)%p
                                for ex,co in rational.terms()},(a,b,c,h),modulus=p)
    ds=[reduce(d.subs(sub)) for d in M['D']]
    w=(2*point[1]+k)%p
    assert int(ds[0].as_expr())%p==w*w%p
    n2=4*w*w*ds[2]-ds[1]**2
    n3=8*w**4*ds[3]-4*w*w*ds[1]*ds[2]+ds[1]**3
    fs=[64*w**6*ds[4]-4*ds[1]*n3-n2**2,
        64*w**8*ds[5]-n2*n3,256*w**10*ds[6]-n3**2]
    assert [f.total_degree() for f in fs]==[4,5,6]
    infinity=[s.Poly(f.as_expr().subs(h,0),a,b,c,modulus=p) for f in fs]
    degree=13;cols=monomials(degree);index={ex:j for j,ex in enumerate(cols)}
    rows=[];labels=[]
    for j,(f,d) in enumerate(zip(infinity,[4,5,6])):
        for multiplier in monomials(degree-d):
            row=[0]*len(cols)
            for ex,co in f.terms():
                target=tuple(x+y for x,y in zip(ex,multiplier))
                row[index[target]]=int(co)%p
            rows.append(row);labels.append([j,list(multiplier)])
    rank,chosen=rank_witness(rows,p)
    assert rank==len(cols)==105
    selected=[rows[i] for i in chosen]
    assert rank_witness(selected,p)[0]==105
    return {'prime':p,'endpoint':{'xi':point[0],'k0':1,'k1':1,'k2':2,'k3':point[1]},
            'w0_mod_p':w,'homogeneous_degrees':[4,5,6],
            'infinity_macaulay_degree':13,'matrix_shape':[len(rows),len(cols)],
            'rank':rank,'selected_independent_rows':[labels[i] for i in chosen],
            'conclusion':'No common projective zero at h=0 for this specialization; hence generic affine degree 120, counting all sheets with multiplicity.'}

if __name__=='__main__':
    M=model();runs=[finite_test(M,p) for p in [101,109]]
    report={'symbolic_checks':{'fixed_target_binding':'PASS','endpoint_double_cover':'PASS',
            'square_root_elimination':'PASS','direction_homogeneity':'PASS','exact_source_deformation':'PASS'},
            'endpoint_equation':str(M['endpoint']),
            'D_coefficients':[str(v) for v in M['D']],
            'B_coefficients':[str(v) for v in M['B']],
            'source_coefficients':[str(v) for v in M['A']],
            'residual_templates':[str(v) for v in M['residual_templates']],
            'runs':runs,
            'local_component_binding':[local_binding(M,p) for p in [101,109]],
            'source_transport':source_transport(),
            'limits':['Presentation localized at w0=2*k3+(5/7)*k0^3 != 0.',
                      'Degree 120 is the total line-incidence degree over the generic endpoint curve; the chosen component has degree at most 120.',
                      'No degree bound over the discriminant coordinate follows from this count alone.',
                      'No global irreducibility, source exactness, rational primitive, or target certificate proved.']}
    (ROOT/'algebraic-curve-model.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'symbolic_checks':report['symbolic_checks'],
          'runs':[{k:r[k] for k in ['prime','endpoint','w0_mod_p','matrix_shape','rank']} for r in runs]},indent=2))
