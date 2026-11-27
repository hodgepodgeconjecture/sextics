"""Generic infinity boundary of the one-component incidence curve.

Exact leading equations and finite-algebra transversality certificates.
Not a global primitive or a HasCertificate proof.
"""
import json
from pathlib import Path
import numpy as np
import sympy as s
from finite_curve_algebra import rref,inverse
ROOT=Path(__file__).resolve().parent

def mon(n):return [(i,d-i) for d in range(n+1) for i in range(d+1)]
def quotient(fs,degree,p):
    cols=mon(degree);idx={e:j for j,e in enumerate(cols)};rows=[]
    for f,d in zip(fs,[6,12]):
        for m in mon(degree-d):
            row=[0]*len(cols)
            for e,co in f.terms():row[idx[tuple(a+b for a,b in zip(e,m))]]=int(co)%p
            rows.append(row)
    rr,piv=rref(rows,p);free=[j for j in range(len(cols)) if j not in piv]
    assert len(free)==72
    nf=np.zeros((72,len(cols)),dtype=np.int64);nf[:,free]=np.eye(72,dtype=np.int64)
    nf[:,piv]=-rr[:len(piv),free].T%p
    return cols,idx,free,nf,len(piv)

def leading():
    eps,t=s.symbols('eps t');b,v,c,h,H,K,L,k=s.symbols('b v c h H K L k')
    x0=K*eps+b*t;x1=eps+v*t;y1=h*t;y2=L*eps+c*t
    B=x0**3+x0*y1*y2+x1*y1
    C=x0**6+x1**6+1+y1**6+y2**6+k*(x0**3*x1*y1+x0*x0*y1*y2*y2+x0*x1*y1**3*y2)
    D=s.Poly(s.expand(k*k*B*B-4*C),t)
    # h=eps*H resolves the vertical plane-sextic component.
    ds=[s.expand(D.nth(j).subs(h,eps*H)) for j in range(7)]
    assert all(ds[j].coeff(eps,i)==0 for j,bound in [(1,5),(2,4),(3,3)] for i in range(bound))
    p4=ds[4].coeff(eps,2);p5=ds[5].coeff(eps,1);p6=ds[6].coeff(eps,0)
    a=s.expand(p4).coeff(H,2);q=s.expand(p4).coeff(H,1);r=s.expand(p4).coeff(H,0)
    u=s.expand(p5).coeff(H,1);z=s.expand(p5).coeff(H,0)
    selector=s.expand(a*z*z-q*u*z+r*u*u)
    assert s.Poly(selector,b,v,c).total_degree()==12
    assert s.expand(p4.subs(H,-z/u)*u*u-selector).cancel()==0
    # Source forcing orders after the blowup. A=(K*eps+b*t)^3,
    # hence A_1,A_2,A_3 have orders 2,1,0. w0_s has order 3.
    source_orders={'W1':2,'W2':1,'W3':0}
    residual_orders=[min(2+3,5+0,4+1),min(1+3,4+0),3+0]
    assert residual_orders==[5,4,3]
    # In the finite h=1 chart, the two kinds of isolated points.
    z0=[s.expand(D.nth(j).subs({eps:0,h:1})) for j in [4,5,6]]
    assert s.expand(z0[0]-k*k*v*v)==0
    assert s.expand(z0[1]-(2*k*(k-2)*b**3*v+2*k*k*b*c*v-4*k*b*b*c*c))==0
    assert s.expand(z0[2]-(k*k*b*b*(b*b+c)**2-4*(b**6+v**6+1+c**6+k*b*v*c)))==0
    # Coefficients forcing the Puiseux balances at the two isolated types.
    first4=s.expand(D.nth(4).subs(h,1)).coeff(eps,1)
    first5=s.expand(D.nth(5).subs(h,1)).coeff(eps,1)
    assert s.expand(first4.subs({v:0,c:0})-2*k*(k-2)*b**3)==0
    assert s.expand(first4.subs(v,0)).coeff(b,1)==-8*K*c*c*k+2*c*k*k
    assert s.expand(first5.subs({b:0,v:0})+24*L*c**5)==0
    # Weighted initial forms at the isolated points, before eliminating
    # the variable with a unit derivative in D6. Corrections from W^2
    # have strictly higher weights (checked below).
    def initial(expr,weights):
        pp=s.Poly(s.expand(expr),eps,b,v,c)
        least=min(sum(x*y for x,y in zip(e,weights)) for e,co in pp.terms())
        return least,s.expand(sum(co*eps**e[0]*b**e[1]*v**e[2]*c**e[3]
            for e,co in pp.terms() if sum(x*y for x,y in zip(e,weights))==least))
    finite=[s.expand(D.nth(j).subs(h,1)) for j in range(7)]
    init=[]
    for weights,expected in [
        ((4,0,2,1),[(4,k*k*v*v+2*k*(k-2)*b**3*eps),
                      (2,2*k*(k-2)*b**3*v-4*k*b*b*c*c)]),
        ((4,2,3,0),[(6,k*k*v*v+2*k*c*(k-4*K*c)*b*eps),
                      (4,-4*k*b*b*c*c-24*L*c**5*eps)])]:
        for j,(order,formula) in zip([4,5],expected):
            got=initial(finite[j],weights)
            assert got[0]==order and s.expand(got[1]-formula)==0
        # w0 is a unit. The recursive W coefficients consequently have
        # at least the following orders, including the W1^2 correction.
        dorders=[initial(finite[j],weights)[0] for j in [1,2,3]]
        w1=dorders[0];w2=min(dorders[1],2*w1);w3=min(dorders[2],w1+w2)
        corrections=[min(w1+w3,2*w2),w2+w3,2*w3]
        assert corrections[0]>expected[0][0] and corrections[1]>expected[1][0]
        # D6 has a simple root in b (first type) or c (second type).
        assert corrections[2]>0
        init.append({'weights_eps_b_v_c':list(weights),
                     'E4_initial':str(expected[0][1]),'E5_initial':str(expected[1][1]),
                     'W_order_lower_bounds':[w1,w2,w3],
                     'residual_correction_order_lower_bounds':corrections})
    # At the specialized generic parameters, no c^6=-1 point has
    # k-4*K*c=0. This excludes the only extra leading-coefficient failure.
    assert s.gcd(s.Poly(c**6+1,c),s.Poly(s.Rational(5,7)-4*c,c)).degree()==0
    return {'variables':(b,v,c,H,K,L,k),'p4':s.expand(p4),'p5':s.expand(p5),
            'p6':s.expand(p6),'U':s.expand(u),'V':s.expand(z),'selector':selector,
            'isolated_weight_checks':init,
            'source_forcing_orders_after_blowup':[3,3,3]}

def certificate(M,p):
    b,v,c,H,K,L,k=M['variables'];subs={K:1,L:2,k:s.Rational(5,7),c:1}
    def poly(expr):
        q=s.Poly(expr,b,v)
        return s.Poly.from_dict({e:int(co.p)*pow(int(co.q),-1,p)%p for e,co in q.terms()},(b,v),modulus=p)
    fs=[poly(M[n].subs(subs)) for n in ['p6','selector']]
    C16,I16,F16,N16,r16=quotient(fs,16,p)
    C17,I17,F17,N17,r17=quotient(fs,17,p)
    basis=[C16[j] for j in F16]
    change=N17[:,[I17[e] for e in basis]];ci=inverse(change,p)
    mats=[]
    for j in range(2):mats.append(ci@N17[:,[I17[tuple(a+int(i==j) for i,a in enumerate(e))] for e in basis]]%p)
    unit=ci@N17[:,I17[(0,0)]]%p;ident=np.eye(72,dtype=np.int64);cache={(0,0):ident}
    def power(e):
        if e not in cache:
            j=next(j for j,x in enumerate(e) if x);a=list(e);a[j]-=1
            cache[e]=mats[j]@power(tuple(a))%p
        return cache[e]
    def op(f):
        out=np.zeros((72,72),dtype=np.int64)
        for e,co in f.terms():out=(out+int(co)*power(e))%p
        return out
    assert np.array_equal(mats[0]@mats[1]%p,mats[1]@mats[0]%p)
    assert all(not np.any(op(f)) for f in fs)
    for j,e in enumerate(basis):assert np.array_equal(power(e)@unit%p,ident[:,j])
    jac=fs[0].diff(b)*fs[1].diff(v)-fs[0].diff(v)*fs[1].diff(b)
    U=poly(M['U'].subs(subs));V=poly(M['V'].subs(subs))
    units={'b':mats[0],'v':mats[1],'selector_Jacobian':op(jac),'U':op(U),'V':op(V)}
    for mat in units.values():assert np.array_equal(mat@inverse(mat,p)%p,ident)
    hmat=-inverse(op(U),p)@op(V)%p
    A=poly(s.expand(M['p4']).coeff(H,2).subs(subs))
    Q=poly(s.expand(M['p4']).coeff(H,1).subs(subs))
    R=poly(s.expand(M['p4']).coeff(H,0).subs(subs))
    assert not np.any((op(A)@hmat@hmat+op(Q)@hmat+op(R))%p)
    # No missing solutions at c=0: resultant of the two binary homogeneous
    # leading parts is nonzero, certified by the degree-17 relation matrix.
    f0=[poly(M[n].subs({K:1,L:2,k:s.Rational(5,7),c:0})) for n in ['p6','selector']]
    rows=[]
    for f,d in zip(f0,[6,12]):
        for i in range(17-d+1):
            row=[0]*18
            for e,co in f.terms():row[e[0]+i]=int(co)%p
            rows.append(row)
    _,piv=rref(rows,p);assert len(piv)==18
    return {'prime':p,'parameters':{'K':1,'L':2,'kappa':'5/7'},
            'quotient_dimension':72,'macaulay_ranks':{'degree16':r16,'degree17':r17},
            'at_c_zero_rank':18,'unit_checks':list(units),
            'checks':{'selector_transverse':'PASS','no_lost_points_at_c_zero':'PASS',
                      'all_vertical_branches_have_b_v_H_nonzero':'PASS','leading_incidence':'PASS'}}

if __name__=='__main__':
    M=leading();results=[certificate(M,p) for p in [1009,1013]]
    report={'leading_equations':{n:str(M[n]) for n in ['p4','p5','p6','U','V']},
      'vertical_source_forcing_orders':M['source_forcing_orders_after_blowup'],
      'isolated_weight_checks':M['isolated_weight_checks'],
      'certificates':results,
      'branch_types':[
        {'name':'vertical-family branches','count':72,'endpoint_ramification':1,'ord_Delta_original':-12,'source_order_lower_bound':3},
        {'name':'c=0, b nonzero','count':6,'endpoint_ramification':4,'ord_Delta_original':-24,'source_order_lower_bound':8},
        {'name':'b=0, c nonzero','count':6,'endpoint_ramification':4,'ord_Delta_original':-20,'Kummer_ramification':3,'source_order_after_Kummer_lower_bound':12}],
      'scope':'Generic endpoint parameters K,L. Formal local balances and source-order bounds are explained in infinity-boundary.md. No finite-boundary classification, compact de Rham reduction, or certificate.'}
    (ROOT/'infinity-boundary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(results,indent=2))
