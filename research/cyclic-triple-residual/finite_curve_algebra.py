"""Construct the 120-dimensional incidence algebra at an exact good fiber.

This verifies a finite-algebra representation, not rational integration.
"""
import json
import time
import math
from pathlib import Path
import numpy as np
import sympy as s
from algebraic_curve_model import model
ROOT=Path(__file__).resolve().parent

def local_residue_checks():
    # For ord(Delta)=m and ramification e=6/gcd(6,m), a term
    # u^n du/r pulls back with exponent e*(n+1)-m/gcd(6,m)-1.
    examples=[]
    for m in range(-12,13):
        g=math.gcd(6,m);e=6//g
        assert ((m//g)%e==0)==(m%6==0)
        if m%6==0:
            n=m//6-1
            assert e*(n+1)-m//g-1==-1
        examples.append({'discriminant_order':m,'ramification_index':e,
                         'residue_allowed_by_character':m%6==0})
    f,fu,fs,xu,xs,xus=s.symbols('f fu fs xu xs xus',nonzero=True)
    old=(fs-fu*xs/xu)*xu
    changed=fs*xu+f*xus-(fu*xs+f*xus)
    assert s.cancel(old-changed)==0
    assert s.cancel(old-(fs*xu-fu*xs))==0
    return {'Kummer_residue_exponent_identity':'PASS',
            'regular_family_coordinate_transport':'PASS','orders_checked':examples,
            'scope':'Local statements on normalized smooth branches. These do not enumerate all actual divisors or settle compact de Rham exactness.'}

def mons(n):
    return [(i,j,k) for d in range(n+1) for i in range(d+1)
            for j in range(d-i+1) for k in [d-i-j]]

def rref(matrix,p):
    a=np.array(matrix,dtype=np.int64)%p;r=0;piv=[]
    for j in range(a.shape[1]):
        found=np.flatnonzero(a[r:,j])
        if not len(found):continue
        q=r+int(found[0]);a[[r,q]]=a[[q,r]]
        a[r]=a[r]*pow(int(a[r,j]),-1,p)%p
        factors=a[:,j].copy();factors[r]=0
        a=(a-factors[:,None]*a[r][None,:])%p
        piv.append(j);r+=1
        if r==a.shape[0]:break
    return a,piv

def inverse(a,p):
    n=len(a);out,piv=rref(np.concatenate([a,np.eye(n,dtype=np.int64)],axis=1),p)
    assert piv[:n]==list(range(n)), 'singular matrix'
    assert np.array_equal(out[:,:n],np.eye(n,dtype=np.int64))
    return out[:,n:]

def quotient(fs,degree,p):
    columns=mons(degree);index={e:i for i,e in enumerate(columns)};rows=[]
    for f,d in zip(fs,[4,5,6]):
        for m in mons(degree-d):
            row=[0]*len(columns)
            for ex,co in f.terms():
                row[index[tuple(x+y for x,y in zip(ex,m))]]=int(co)%p
            rows.append(row)
    rr,piv=rref(rows,p);free=[j for j in range(len(columns)) if j not in piv]
    assert len(free)==120
    normal=np.zeros((120,len(columns)),dtype=np.int64)
    normal[:,free]=np.eye(120,dtype=np.int64)
    normal[:,piv]=-rr[:len(piv),free].T%p
    return columns,index,free,normal,len(piv)

def run(p=1009):
    start=time.monotonic();M=model();a,b,c,h=M['vars'];xi,k0,k2,k3=M['params']
    kap=5*pow(7,-1,p)%p
    point=next((x,z) for x in range(1,p) for z in range(p)
        if (x**6+1+1+2**6+z*z+kap*z)%p==0 and (2*z+kap)%p)
    sub={xi:point[0],k0:1,k2:2,k3:point[1],h:1}
    variables=(a,b,c)
    def poly(expr):
        q=s.Poly(expr,variables)
        return s.Poly.from_dict({e:int(v.p)*pow(int(v.q),-1,p)%p for e,v in q.terms()},variables,modulus=p)
    ds=[poly(q.subs(sub)) for q in M['D']];bs=[poly(q.subs(sub)) for q in M['B']]
    w=(2*point[1]+kap)%p
    assert int(ds[0].as_expr())%p==w*w%p
    n2=4*w*w*ds[2]-ds[1]**2;n3=8*w**4*ds[3]-4*w*w*ds[1]*ds[2]+ds[1]**3
    fs=[64*w**6*ds[4]-4*ds[1]*n3-n2**2,64*w**8*ds[5]-n2*n3,256*w**10*ds[6]-n3**2]
    C12,I12,F12,N12,rank12=quotient(fs,12,p)
    print(f'p={p} degree12 rank={rank12}, quotient=120',flush=True)
    C13,I13,F13,N13,rank13=quotient(fs,13,p)
    basis=[C12[j] for j in F12]
    change=N13[:,[I13[e] for e in basis]];ci=inverse(change,p)
    mats=[]
    for j in range(3):
        shifted=[tuple(v+int(i==j) for i,v in enumerate(e)) for e in basis]
        mats.append(ci@N13[:,[I13[e] for e in shifted]]%p)
    unit=ci@N13[:,I13[(0,0,0)]]%p
    ident=np.eye(120,dtype=np.int64)
    assert all(np.array_equal(x@y%p,y@x%p) for x in mats for y in mats)
    cache={(0,0,0):ident}
    def monomial(e):
        if e not in cache:
            j=next(j for j,v in enumerate(e) if v)
            prev=list(e);prev[j]-=1
            cache[e]=mats[j]@monomial(tuple(prev))%p
        return cache[e]
    def operator(f):
        ans=np.zeros((120,120),dtype=np.int64)
        for e,co in f.terms():ans=(ans+(int(co)%p)*monomial(e))%p
        return ans
    assert all(not np.any(operator(f)) for f in fs)
    # Every advertised basis monomial applied to 1 returns its basis vector.
    for j,e in enumerate(basis):assert np.array_equal(monomial(e)@unit%p,ident[:,j])
    Z=[poly(point[1]),(ds[1]*pow(2*w,-1,p)-kap*bs[1])*pow(2,-1,p),
       (n2*pow(8*w**3,-1,p)-kap*bs[2])*pow(2,-1,p),
       (n3*pow(16*w**5,-1,p)-kap*bs[3])*pow(2,-1,p)]
    z0,z1,z2,z3=Z
    Delta=z2**2*z1**2-4*z3*z1**3-4*z2**3*z0-27*z3**2*z0**2+18*z3*z2*z1*z0
    Md=operator(Delta)
    print(f'p={p} multiplication matrices verified; constructing discriminant powers',flush=True)
    vectors=[];v=unit
    for _ in range(120):vectors.append(v);v=Md@v%p
    krylov=np.column_stack(vectors);ki=inverse(krylov,p)
    coefficients=ki@v%p
    assert np.array_equal(krylov@coefficients%p,v)
    T=s.Symbol('T');minimal=s.Poly(T**120-sum(int(co)*T**j for j,co in enumerate(coefficients)),T,modulus=p)
    assert s.gcd(minimal,minimal.diff()).degree()==0
    assert minimal.nth(0)!=0
    # All slopes can now be represented by polynomials in Delta on this fiber.
    slope_repr=[(ki@(m@unit%p)%p).tolist() for m in mats]
    for m,co in zip(mats,slope_repr):assert np.array_equal(krylov@np.array(co,dtype=np.int64)%p,m@unit%p)
    # Differentiate the reduced incidence in xi and in the actual source
    # deformation. Solve in the full finite algebra, without choosing roots.
    wx=-12*point[0]**5*pow(w,-1,p)%p
    ws=2*point[0]**3*pow(w,-1,p)%p  # k0=1 on this test fiber
    Dxi=[poly(s.diff(q,xi).subs(sub)) for q in M['D']]
    Bxi=[poly(s.diff(q,xi).subs(sub)) for q in M['B']]
    As=[poly(q.subs(sub)) for q in M['A']]
    tv=M['template_vars'];template_sub=dict(zip(tv,[w]+[q.as_expr() for q in ds[1:]]))
    forcing=[[],[]]
    for template in M['residual_templates']:
        fw=poly(template.diff(tv[0]).subs(template_sub,simultaneous=True))
        fd=[poly(template.diff(q).subs(template_sub,simultaneous=True)) for q in tv[1:]]
        fx=fw*wx;fsrc=fw*ws
        for j,q in enumerate(fd,1):fx+=q*Dxi[j];fsrc+=q*4*As[j]
        forcing[0].append(operator(fx)@unit%p)
        forcing[1].append(operator(fsrc)@unit%p)
    jac=np.block([[operator(f.diff(x)) for x in variables] for f in fs])
    rhs=-np.column_stack([np.concatenate(f) for f in forcing])%p
    sol,piv=rref(np.concatenate([jac,rhs],axis=1),p)
    assert piv[:360]==list(range(360))
    velocities=sol[:,-2:]
    assert np.array_equal(jac@velocities%p,rhs)
    def delta_explicit(dD,dB,dw):
        dn2=8*w*dw*ds[2]+4*w*w*dD[2]-2*ds[1]*dD[1]
        dn3=(32*w**3*dw*ds[3]+8*w**4*dD[3]-8*w*dw*ds[1]*ds[2]
              -4*w*w*(dD[1]*ds[2]+ds[1]*dD[2])+3*ds[1]**2*dD[1])
        dW=[poly(dw),dD[1]*pow(2*w,-1,p)-ds[1]*dw*pow(2*w*w,-1,p),
            dn2*pow(8*w**3,-1,p)-3*n2*dw*pow(8*w**4,-1,p),
            dn3*pow(16*w**5,-1,p)-5*n3*dw*pow(16*w**6,-1,p)]
        dZ=[(q-kap*b)*pow(2,-1,p) for q,b in zip(dW,dB)]
        zz=s.symbols('z0:4')
        disc=zz[2]**2*zz[1]**2-4*zz[3]*zz[1]**3-4*zz[2]**3*zz[0]-27*zz[3]**2*zz[0]**2+18*zz[3]*zz[2]*zz[1]*zz[0]
        out=poly(0)
        for j in range(4):
            grad=poly(s.diff(disc,zz[j]).subs(dict(zip(zz,[z.as_expr() for z in Z])),simultaneous=True))
            out+=grad*dZ[j]
        return out
    explicit=[delta_explicit(Dxi,Bxi,wx),delta_explicit([4*q for q in As],[poly(0)]*4,ws)]
    delta_variations=[]
    for col in range(2):
        v=operator(explicit[col])@unit%p
        for j,x in enumerate(variables):v=(v+operator(Delta.diff(x))@velocities[j*120:(j+1)*120,col])%p
        delta_variations.append(v)
    delta_xi,delta_dot=delta_variations
    mdi=inverse(Md,p)
    source=(2*point[0]**3*pow(3*w**3,-1,p)*unit+
             pow(18*w,-1,p)*(mdi@delta_dot%p))%p
    old_ds=(delta_dot+pow(6*point[0]**2,-1,p)*delta_xi)%p
    old_density=(pow(9*point[0]**3*w,-1,p)*unit+pow(18*w,-1,p)*(mdi@old_ds%p))%p
    correction=((-pow(9*point[0]**3*w,-1,p)+2*point[0]**3*pow(3*w**3,-1,p))*unit
                -pow(108*point[0]**2*w,-1,p)*(mdi@delta_xi%p))%p
    assert np.array_equal((old_density+correction)%p,source)
    # Derivative of each Delta^j/r is (j-1/6) Delta^(j-1) Delta_xi/r.
    connection_columns=[];v=mdi@delta_xi%p
    for j in range(120):
        connection_columns.append(ki@((j-pow(6,-1,p))*v%p)%p)
        v=Md@v%p
    connection=np.column_stack(connection_columns)
    # Product rule for multiplication by Delta, including reduction of its
    # top power via the minimal polynomial.
    companion=ki@Md@krylov%p
    powers=[];v=delta_xi
    for j in range(120):powers.append(v);v=Md@v%p
    dp=120*powers[119]%p
    for j in range(1,120):dp=(dp-j*int(coefficients[j])*powers[j-1])%p
    coefficient_derivatives=ki@dp%p
    dcompanion=np.zeros((120,120),dtype=np.int64);dcompanion[:,-1]=coefficient_derivatives
    mult_delta_xi=ki@np.column_stack(powers)%p
    assert np.array_equal((dcompanion+connection@companion-companion@connection)%p,mult_delta_xi)
    print(f'p={p} full source and differential-field connection reconstructed',flush=True)
    result={'prime':p,'endpoint':{'xi':point[0],'k0':1,'k1':1,'k2':2,'k3':point[1]},
      'w0_mod_p':w,'macaulay_ranks':{'degree_12':rank12,'degree_13':rank13},
      'quotient_dimension':120,'basis_exponents':basis,
      'discriminant_minimal_polynomial_coefficients_ascending':[int(minimal.nth(j))%p for j in range(121)],
      'slope_polynomials_in_discriminant':dict(zip(['a','b','c'],slope_repr)),
      'discriminant_xi_derivative_coefficients':(ki@delta_xi%p).tolist(),
      'discriminant_source_derivative_coefficients':(ki@delta_dot%p).tolist(),
      'transported_source_coefficients':(ki@source%p).tolist(),
      'reduction_coefficients_xi_derivative':coefficient_derivatives.tolist(),
      'connection_reconstruction':'column j = (j-1/6) times Delta^(j-1)*Delta_xi, reduced modulo the reported polynomial; j=0 uses Delta inverse',
      'checks':{'multiplication_commutes':'PASS','incidence_relations':'PASS','monomial_basis_binding':'PASS',
                'discriminant_is_cyclic_generator':'PASS','minimal_polynomial_squarefree':'PASS',
                'discriminant_invertible':'PASS','slope_reconstruction':'PASS',
                'implicit_parameter_and_source_derivatives':'PASS','actual_source_transport':'PASS',
                'connection_product_rule':'PASS'},
      'elapsed_seconds':round(time.monotonic()-start,3),
      'scope':'Exact specialized finite algebra. Does not give global parameter-dependent polynomial coefficients, ramification divisors, residues, rational primitive, or HasCertificate.'}
    return result

if __name__=='__main__':
    results=[run(p) for p in [1009,1013]]
    (ROOT/'finite-curve-algebra.json').write_text(json.dumps({'fibers':results,
                       'local_residue_checks':local_residue_checks()},indent=2)+'\n')
    print(json.dumps([{k:r[k] for k in ['prime','quotient_dimension','macaulay_ranks','checks','elapsed_seconds']} for r in results],indent=2))
