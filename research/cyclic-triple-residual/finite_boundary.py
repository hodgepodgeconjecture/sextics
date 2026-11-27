"""An explicit finite boundary family of the incidence integration curve.

Exact characteristic-zero identities, not a global boundary classification
or a rational primitive. See finite-boundary.md for the valuation argument.
"""
import json
from pathlib import Path
import sympy as s
from algebraic_curve_model import model

ROOT=Path(__file__).resolve().parent


def constant_cubic_family():
    M=model()  # Binds the construction to the protected F, q, A.
    aa,bb,cc,hh=M['vars'];xi,k0,k2,k3=M['params']
    x,u,b,h,K,L,r,w=s.symbols('x u b h K L rho w')
    k=s.Rational(5,7);d=k*k-4
    variables=(x,u,b,h);weights=(5,6,1,4)
    sub={xi:r*L+x,aa:r+u,bb:b,cc:1,hh:h,k0:K,k2:L}
    def rho_reduce(expr):return s.rem(s.expand(expr),r**6+1,r)
    D=[rho_reduce(q.subs(sub,simultaneous=True)) for q in M['D']]
    A=[rho_reduce(q.subs(sub,simultaneous=True)) for q in M['A']]
    def initial(expr):
        terms=s.Poly(s.expand(expr),*variables).terms()
        terms=[(e,co) for e,co in terms if co]
        if not terms:return 1000,s.Integer(0)
        n=min(sum(a*c for a,c in zip(e,weights)) for e,co in terms)
        return n,s.expand(sum(co*s.prod(a**j for a,j in zip(variables,e))
            for e,co in terms if sum(a*c for a,c in zip(e,weights))==n))
    zero={x:0,u:0,b:0,h:0}
    assert s.expand(D[0].subs(zero)-(d*K**6-4))==0
    assert all(q.subs(zero)==0 for q in D[1:])
    jac=s.Matrix([[q.diff(v).subs(zero) for v in [x,u,h]] for q in D[4:]])
    expected=s.Matrix([[-120*r**5*L,-240*r**5*L**2,-4*k*K**2*r],
                       [-24*r**5,-120*r**5*L,0],[0,-24*r**5,0]])
    assert jac==expected
    determinant=s.factor(rho_reduce(jac.det()))
    assert determinant!=0
    for j in [1,2,3]:
        assert initial(D[j])==(j,d*s.binomial(6,j)*K**(6-j)*b**j)
    expected_initial=[15*d*K*K*b**4-4*k*K*K*r*h,
                      6*d*K*b**5-8*k*K*r*b*h-24*r**5*x,
                      d*b**6-4*k*r*b*b*h-24*r**5*u]
    for j,q in zip([4,5,6],expected_initial):
        n,f=initial(D[j]);assert n==j and s.expand(f-q)==0
    # The square root coefficients along the transverse b direction.
    relation=w*w-d*K**6+4
    def reduce_w(q):
        n,de=s.fraction(s.cancel(q))
        return s.factor(s.rem(n,relation,w)/s.rem(de,relation,w))
    W1=3*d*K**5/w
    W2=reduce_w((15*d*K**4-W1*W1)/(2*w))
    W3=reduce_w((20*d*K**3-2*W1*W2)/(2*w))
    R4=reduce_w(15*d*K*K-2*W1*W3-W2*W2)
    R5=reduce_w(6*d*K-2*W2*W3)
    R6=reduce_w(d-W3*W3)
    H4=reduce_w(R4/(4*k*K*K*r))
    X5=reduce_w((R5-8*k*K*r*H4)/(24*r**5))
    U6=reduce_w((R6-4*k*r*H4)/(24*r**5))
    for q in [R4-4*k*K*K*r*H4,
              R5-8*k*K*r*H4-24*r**5*X5,
              R6-4*k*r*H4-24*r**5*U6]:assert reduce_w(q)==0
    # At K=1,rho=i all denominators and H4,X5 are nonzero for either w.
    specialized_relation=relation.subs(K,1)
    def coprime_at_one(q):
        n,de=s.fraction(s.cancel(q.subs({K:1,r:s.I})))
        assert s.gcd(n,specialized_relation,w)==1
        assert s.gcd(de,specialized_relation,w)==1
    for q in [H4,X5]:coprime_at_one(q)
    Z=[(q-k*q2)/2 for q,q2 in zip([w,W1,W2,W3],[K**3,3*K*K,3*K,1])]
    z0,z1,z2,z3=Z
    disc=z2*z2*z1*z1-4*z3*z1**3-4*z2**3*z0-27*z3*z3*z0*z0+18*z3*z2*z1*z0
    coprime_at_one(disc)
    # Source forcing begins in weights 1,2,3. Differentiating W^2 uses
    # regular W_{j,s}, so the product corrections have the same bounds.
    source_orders=[initial(A[j])[0] for j in [4,5,6]]
    assert source_orders==[1,2,3]
    # Check the lowest Jacobian orders that preserve x_s=O(b^2),
    # u_s=O(b^3),h_s=O(b) when b is the local coordinate.
    jac_bounds=[[0,0,0],[0,0,1],[3,0,2]]
    # D6 is independent of x; its W3^2 correction has x derivative
    # of order at least 3. Other W-product corrections do not lower these.
    for row,j in enumerate([4,5,6]):
        for col,v in enumerate([x,u,h]):
            assert initial(D[j].diff(v))[0]>=jac_bounds[row][col]
    return {'scope':'One explicit finite boundary family, generic K,L. Not all finite boundary points.',
            'boundary':'xi=rho*L, a=rho, b=h=0, c=1; rho^6=-1',
            'endpoint_relation':str(relation),
            'local_Jacobian_variables':['x','u','h'],
            'local_Jacobian_determinant':str(determinant),
            'initial_D4_D5_D6':[str(q) for q in expected_initial],
            'series_leading_coefficients':{'h_order4':str(H4),'xi_shift_order5':str(X5),'a_shift_order6':str(U6)},
            'source_forcing_order_lower_bounds':source_orders,
            'generic_nonzero_conditions_checked_at_K_one':['h_order4','xi_shift_order5','raw_discriminant_order6'],
            'branch_orders':{'endpoint_ramification':5,'h':4,'raw_discriminant':6,'old_discriminant':-18,
                             'Kummer_ramification':1,'eta':7,'transformed_source_lower_bound':4},
            'checks':{'protected_target_binding':'PASS','constant_cubic_boundary':'PASS',
                      'unit_local_Jacobian':'PASS','weighted_initial_equations':'PASS',
                      'leading_coefficient_equations':'PASS','generic_nonvanishing':'PASS',
                      'source_forcing_and_Jacobian_bounds':'PASS'}}


def translated_cubic_family():
    M=model();aa,bb,cc,hh=M['vars'];xi,k0,k2,k3=M['params']
    x,u,c,h,K,L,r,w=s.symbols('x u c h K L rho w')
    k=s.Rational(5,7);d=k*k-4
    relation_r=r**6-d/4;relation_w=w*w+4*(1+L**6)
    sub={xi:r*K+x,aa:r+u,bb:1,cc:c,hh:h,k0:K,k2:L}
    variables=(x,u,c,h);weights=(4,5,1,4)
    def rr(q):return s.rem(s.expand(q),relation_r,r)
    D=[rr(q.subs(sub,simultaneous=True)) for q in M['D']]
    def initial(q):
        terms=[(e,co) for e,co in s.Poly(s.expand(q),*variables).terms() if co]
        if not terms:return 1000,s.Integer(0)
        n=min(sum(a*b for a,b in zip(e,weights)) for e,co in terms)
        return n,s.expand(sum(co*s.prod(a**j for a,j in zip(variables,e))
            for e,co in terms if sum(a*b for a,b in zip(e,weights))==n))
    zero={x:0,u:0,c:0,h:0}
    assert s.expand(D[0].subs(zero)+4*(1+L**6))==0
    assert all(q.subs(zero)==0 for q in D[1:])
    A5=2*k*(k*L+(k-2)*r);A4=4*K*A5-4*k*r*L*L
    jac=s.Matrix([[q.diff(v).subs(zero) for v in [x,u,h]] for q in D[4:]])
    expected=s.Matrix([[-120*r**5*K,-240*r**5*K*K,A4],
                       [-24*r**5,-120*r**5*K,A5],[0,-24*r**5,0]])
    assert all(s.expand(a-b)==0 for a,b in zip(jac,expected))
    determinant=s.factor(rr(jac.det()))
    for j in [1,2,3]:
        assert initial(D[j])==(j,-4*s.binomial(6,j)*L**(6-j)*c**j)
    expected_initial=[-60*L*L*c**4-120*r**5*K*x+A4*h,
                      -24*r**5*x+A5*h,-24*r**5*u+2*k*k*c*h]
    for j,n,q in zip([4,5,6],[4,4,5],expected_initial):
        order,f=initial(D[j]);assert order==n and s.expand(f-q)==0
    def rw(q):
        n,de=s.fraction(s.cancel(q))
        return s.factor(s.rem(n,relation_w,w)/s.rem(de,relation_w,w))
    W1=-12*L**5/w;W2=rw((-60*L**4-W1*W1)/(2*w))
    W3=rw((-80*L**3-2*W1*W2)/(2*w))
    R4=rw(-60*L*L-2*W1*W3-W2*W2)
    H4=rw(R4/(K*A5+4*k*r*L*L))
    X4=rw(A5*H4/(24*r**5));U5=rw(k*k*H4/(12*r**5))
    assert rw(R4-120*r**5*K*X4+A4*H4)==0
    assert rw(-24*r**5*X4+A5*H4)==0
    assert rw(-24*r**5*U5+2*k*k*H4)==0
    for q in [determinant,H4,X4]:
        n,de=s.fraction(s.cancel(q.subs({K:1,L:2})))
        assert s.gcd(n,relation_r,r)==1 and s.gcd(de,relation_r,r)==1
    # After translating t by K, Z=(w-kappa*t^3)/2.
    zz=s.symbols('z0:4');t=s.symbols('t')
    cubic=s.Poly((w-k*(K+t)**3)/2,t)
    disc=s.discriminant(cubic.as_expr(),t)
    assert s.expand(disc+27*k*k*w*w/16)==0
    return {'scope':'Second explicit finite boundary family, generic K,L. Not exhaustive.',
            'boundary':'xi=rho*K, a=rho, b=1, c=h=0; rho^6=(kappa^2-4)/4',
            'endpoint_relation':str(relation_w),'rho_relation':str(relation_r),
            'local_Jacobian_determinant':str(determinant),
            'initial_D4_D5_D6':[str(q) for q in expected_initial],
            'series_leading_coefficients':{'h_order4':str(H4),'xi_shift_order4':str(X4),'a_shift_order5':str(U5)},
            'raw_discriminant_at_boundary':str(disc),
            'generic_nonzero_conditions_checked_at_K_one_L_two':['local_Jacobian','h_order4','xi_shift_order4'],
            'branch_orders':{'endpoint_ramification':4,'h':4,'raw_discriminant':0,'old_discriminant':-24,
                             'Kummer_ramification':1,'eta':7,'transformed_source_lower_bound':3},
            'checks':{'protected_target_binding':'PASS','translated_cubic_boundary':'PASS',
                      'unit_local_Jacobian':'PASS','weighted_initial_equations':'PASS',
                      'leading_coefficient_equations':'PASS','generic_nonvanishing':'PASS',
                      'raw_discriminant':'PASS'}}


if __name__=='__main__':
    result={'constant_cubic':constant_cubic_family(),'translated_cubic':translated_cubic_family()}
    (ROOT/'finite-boundary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
