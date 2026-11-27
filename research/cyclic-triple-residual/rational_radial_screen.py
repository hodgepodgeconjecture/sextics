"""Screen a pole-informed radial rational-current span on full incidence points."""
import json
import random
from pathlib import Path
import sympy as s
from fractions import Fraction
from source_root_chart import Adic,build
from line_exactness import evaluate,inverse,mv,compile_poly,derivative,mul,power
ROOT=Path(__file__).resolve().parent

def columns():
    # u=k0/k1, v=k2/k1, w=k3/k1^3. These are precisely the six
    # monomial types of the characteristic-zero leading radial current.
    bases=[(3,0,0),(2,2,-1),(1,3,-1),(1,2,0),(0,3,0),(-1,3,1)]
    factors=[None,0,1,2,3,4,5,6]
    weights=[1,0,0,0,2,1,0,1,1,1,3]
    out=[]
    for u,v,w in bases:
        for factor in factors:
            e=[0]*11;e[7]=u;e[9]=v;e[10]=w;e[8]=-u-v-3*w
            if factor is not None:
                e[factor]+=1;e[8]-=weights[factor]
            assert sum(a*b for a,b in zip(e,weights))==0
            out.append(tuple(e))
    assert len(set(out))==48
    return out

def row(p,n,direction,data,cols,independent=False):
    Q,D,dc,dg=data;mod=p**n
    ri=next(i for i in range(p) if (i*i+1)%p==0)
    point=[ri,0,0,ri,0,0,0,p*direction[0],1,p*direction[1],p*direction[2]]
    for _ in range(n.bit_length()+2):
        jac=[[evaluate(v,point,mod) for v in r[:7]] for r in D]
        correction=mv(inverse(jac,mod),[evaluate(q,point,mod) for q in Q],mod)
        point[:7]=[(a-b)%mod for a,b in zip(point[:7],correction)]
    assert all(evaluate(q,point,mod)==0 for q in Q)
    jac=[[evaluate(v,point,mod) for v in r] for r in D]
    ji=inverse([r[:7] for r in jac],mod)
    # The unweighted radial field on k1=1, not the affine-cone Euler field.
    Rdep=mv(ji,[-sum(r[j]*point[j] for j in [7,9,10]) for r in jac],mod)
    R=[Adic(v,n) for v in Rdep+[point[7],0,point[9],point[10]]]
    fields=[R]
    if independent:
        fields=[]
        for j in [7,9,10]:
            v=mv(ji,[-r[j]*point[j] for r in jac],mod)
            fields.append([Adic(q,n) for q in v+[point[l] if l==j else 0 for l in range(7,11)]])
    A=mul(power([point[7],point[2]],3,mod),power([point[0],point[1]],3,mod),mod)
    A+=[0]*(7-len(A))
    U=mv(ji,[7*v for v in A],mod)
    Z=[point[6],point[5],point[4],point[10]]
    Delta=Adic(evaluate(dc,Z,mod),n)
    grad=[Adic(evaluate(v,Z,mod),n) for v in dg]
    Ds=sum(a*Adic(b,n) for a,b in zip(grad,[U[6],U[5],U[4],0]))
    point=[Adic(v,n) for v in point];xi=point[0]
    def monomial(exps):
        out=Adic(1)
        for x,e in zip(point,exps):
            if e:out=out*x**e
        return out
    vals=[]
    for R in fields:
        RD=sum(a*b for a,b in zip(grad,[R[6],R[5],R[4],R[10]]))
        div=1 if independent else 3
        for exps in cols:
            g=monomial(exps);Rg=Adic(0)
            for j,e in enumerate(exps):
                if e:
                    reduced=list(exps);reduced[j]-=1
                    Rg+=e*monomial(reduced)*R[j]
            vals.append(18*Delta*xi*Rg+(18*div*Delta*xi-90*Delta*R[0]-3*xi*RD)*g)
    rhs=-2*Delta*point[7]**3/xi**2-xi*Ds
    return vals+[rhs]

def verify_formula():
    D,RD,x,Rx,g,Rg,Ds,k,nu=s.symbols('D RD x Rx g Rg Ds k nu')
    J=1/(6*x**5);RJ=s.diff(J,x)*Rx
    lhs=RJ*g+J*Rg+nu*J*g-J*g*RD/(6*D)
    assert s.cancel(18*D*x*lhs/J-(18*D*x*Rg+(18*nu*D*x-90*D*Rx-3*x*RD)*g))==0
    B=-k**3/(54*x**8)-J*Ds/(18*D)
    assert s.cancel(18*D*x*B/J-(-2*D*k**3/x**2-x*Ds))==0
    alpha,beta,gamma,kappa=s.symbols('alpha beta gamma k')
    exact=json.loads((ROOT/'source-deformation-exact.json').read_text())
    C=s.sympify(exact['compact_C0_lambda3'],locals={'alpha':alpha,'beta':beta,'gamma':gamma,'k':kappa})
    monomials=[alpha**3,alpha**2*beta**2/gamma,alpha*beta**3/gamma,
               alpha*beta**2,beta**3,beta**3*gamma/alpha]
    coeffs=[-5*s.I/288,3/(4*kappa**3),3*s.I/(4*kappa**3),
            3*s.I/kappa**4,-12/kappa**4,9*s.I/(8*kappa**3)]
    assert s.cancel(sum(a*b for a,b in zip(coeffs,monomials))-s.I*C/8)==0

def screen(p,n=48,independent=False):
    Adic.p=p;Adic.cap=2*n
    Q,D,_=build('X1');t=s.Symbol('t');z=s.symbols('a b c d')
    dc=compile_poly(s.discriminant(z[0]*t**3+z[1]*t*t+z[2]*t+z[3],t),z)
    data=(Q,D,dc,[derivative(dc,j) for j in range(4)])
    cols=columns();count=len(cols)*(3 if independent else 1)
    basis=[];unused=set(range(count));rng=random.Random(20261127)
    samples=[];pivots=[];witness=None
    for sample in range(count+20):
        direction=[rng.randrange(1,p) for _ in range(3)];samples.append(direction)
        r=row(p,n,direction,data,cols,independent)
        for j,b in basis:
            factor=r[j];r=[x-factor*y for x,y in zip(r,b)];r[j]=Adic(0)
        choices=[j for j in unused if r[j]!=0]
        if choices:
            j=min(choices,key=lambda j:r[j].val())
            assert min(r[l].n for l in unused)>r[j].val()
            valuation=r[j].val();iv=r[j].inv();r=[x*iv for x in r];r[j]=Adic(1)
            basis.append((j,r));unused.remove(j)
            pivots.append({'sample':sample,'column':j,'valuation':valuation,'precision':min(x.n for x in r)})
        elif not unused and r[-1]!=0:
            v=r[-1].val()
            witness={'sample':sample,'valuation':v,'absolute_precision':r[-1].n,
                     'unit_mod_p':int(r[-1].q/Fraction(p)**v)%p}
            break
    return {'prime':p,'input_precision':n,'independent_components':independent,
            'columns':count,'rank_lower_bound':len(basis),
            'directions':samples,'pivots':pivots,'inconsistency_witness':witness,
            'status':'EXCLUDED' if witness else 'INCONCLUSIVE',
            'scope':'Only b_i=J*x_i*g_i on k1=1, with g_i in the stated 48-dimensional Laurent-monomial span; g_i are equal in radial mode. No exclusion of other rational currents.'}

if __name__=='__main__':
    verify_formula();runs=[screen(p,independent=ind) for ind in [False,True] for p in [101,109]]
    report={'cleared_equation_identity':'PASS','column_exponents':columns(),'runs':runs}
    (ROOT/'rational-radial-screen.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps([{k:r[k] for k in ['prime','columns','rank_lower_bound','status','inconsistency_witness']} for r in runs],indent=2))
