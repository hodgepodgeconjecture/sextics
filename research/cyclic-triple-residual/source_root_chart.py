"""Move weighted-line infinity to X0=0 and rebuild the actual source.

Precision-tracked p-adic arithmetic retains poles in the coefficient
Jacobian inverse. This is an exact-congruence experiment, not a Lean proof.
The second chart has X0=k0 constant and X1=k1+b*t.
"""
import argparse
import json
from fractions import Fraction
from pathlib import Path
import sympy as sp
from line_exactness import (compile_poly,derivative,evaluate,inverse,mv,
                            add,mul,scale,power,diff,trim,det_poly)


class Adic:
    p=101
    cap=60
    def __init__(self,x=0,n=None):
        if isinstance(x,Adic):
            self.q=x.q;self.n=min(x.n,n) if n is not None else x.n;return
        self.n=self.cap if n is None else n
        q=Fraction(x)
        if not q: self.q=q;return
        num,den=q.numerator,q.denominator;v=0
        while num % self.p==0: num//=self.p;v+=1
        while den % self.p==0: den//=self.p;v-=1
        if v>=self.n:self.q=Fraction(0);return
        modulus=self.p**(self.n-v)
        unit=num*pow(den,-1,modulus)%modulus
        self.q=Fraction(unit*self.p**v) if v>=0 else Fraction(unit,self.p**(-v))
    def val(self):
        if not self.q:return self.n
        num,den=self.q.numerator,self.q.denominator;v=0
        while num % self.p==0:num//=self.p;v+=1
        while den % self.p==0:den//=self.p;v-=1
        return v
    def coerce(self,o):return o if isinstance(o,Adic) else Adic(o)
    def __add__(self,o):
        o=self.coerce(o);return Adic(self.q+o.q,min(self.n,o.n))
    __radd__=__add__
    def __neg__(self):return Adic(-self.q,self.n)
    def __sub__(self,o):return self+-self.coerce(o)
    def __rsub__(self,o):return self.coerce(o)+-self
    def __mul__(self,o):
        o=self.coerce(o)
        return Adic(self.q*o.q,min(self.n+o.val(),o.n+self.val()))
    __rmul__=__mul__
    def inv(self):
        assert self.q,'division by unresolved zero'
        return Adic(1/self.q,self.n-2*self.val())
    def __truediv__(self,o):return self*self.coerce(o).inv()
    def __rtruediv__(self,o):return self.coerce(o)*self.inv()
    def __pow__(self,n,modulus=None):
        if n<0:return self.inv()**(-n)
        out=Adic(1)
        for _ in range(n):out=out*self
        return out
    def __mod__(self,modulus):return self
    def __eq__(self,o):return self.q==self.coerce(o).q
    def __repr__(self):return f'{self.q}+O({self.p}^{self.n})'


def inv_adic(a):
    size=len(a)
    m=[[Adic(v) for v in row]+[Adic(int(i==j)) for j in range(size)] for i,row in enumerate(a)]
    for j in range(size):
        pivot=min((i for i in range(j,size) if m[i][j]!=0),key=lambda i:m[i][j].val())
        m[j],m[pivot]=m[pivot],m[j]
        u=m[j][j].inv();m[j]=[v*u for v in m[j]]
        for i in range(size):
            if i!=j:
                v=m[i][j];m[i]=[x-v*y for x,y in zip(m[i],m[j])]
    return [row[size:] for row in m]


def build(anchor):
    t=sp.Symbol('t')
    xi,a,b,c,d,e,f,k0,k1,k2,k3=vs=sp.symbols('xi a b c d e f k0 k1 k2 k3')
    y0=xi+a*t;y1=t;y2=k2+c*t;z=k3+d*t+e*t*t+f*t**3
    x0=k0+b*t if anchor=='X1' else k0
    x1=k1 if anchor=='X1' else k1+b*t
    H7=7*(y0**6+y1**6+x0**6+x1**6+y2**6+z*z)+5*(x0**3*x1*y0*y1+x0**3*z+
        x0*x0*y0*y1*y2*y2+x0*x1*y1**3*y2+x0*z*y1*y2+x1*z*y0*y1)
    poly=sp.Poly(sp.expand(H7),t)
    Q=[compile_poly(poly.nth(i),vs) for i in range(7)]
    D=[[derivative(q,j) for j in range(11)] for q in Q]
    DD=[[[derivative(row[j],l) for l in range(11)] for j in range(11)] for row in D]
    return Q,D,DD


def main(p,n,output):
    Adic.p=p;Adic.cap=2*n
    modulus=p**n
    Q,D,_=build('X1')
    ri=next(i for i in range(p) if (i*i+1)%p==0)
    point=[ri,0,0,ri,0,0,0,p,1,2*p,3*p]
    for _ in range(n.bit_length()+2):
        jac=[[evaluate(d,point,modulus) for d in row[:7]] for row in D]
        delta=mv(inverse(jac,modulus),[evaluate(q,point,modulus) for q in Q],modulus)
        point[:7]=[(x-v)%modulus for x,v in zip(point[:7],delta)]
    assert all(evaluate(q,point,modulus)==0 for q in Q)
    point=[Adic(v,n) for v in point]
    xi,a,b,c,d,e,f,k0,k1,k2,k3=point
    h=b/k0
    transformed=[xi,a-xi*h,-k1*h,c-k2*h,d-3*k3*h,e-2*d*h+3*k3*h*h,
                 f-e*h+d*h*h-k3*h**3,k0,k1,k2,k3]
    Q,D,DD=build('X0')
    vals=[Adic(evaluate(q,transformed,1)) for q in Q]
    assert all(v==0 for v in vals)
    jac=[[Adic(evaluate(d,transformed,1)) for d in row] for row in D]
    ji=inv_adic([row[:7] for row in jac])
    identity=[[sum(jac[i][k]*ji[k][j] for k in range(7))-int(i==j) for j in range(7)] for i in range(7)]
    assert all(v==0 for row in identity for v in row)
    eps=mv(ji,[7,0,0,0,0,0,0],1)+[Adic(0)]*4
    directions=[eps]
    for j in [8,9,10]:
        dy=mv(ji,[-row[j] for row in jac],1)
        directions.append(dy+[Adic(int(j==l)) for l in range(7,11)])
    hess=[[[Adic(evaluate(dd,transformed,1)) for dd in row] for row in qq] for qq in DD]
    second=[]
    for v in directions:
        forcing=[sum(h[i][j]*eps[i]*v[j] for i in range(11) for j in range(11)) for h in hess]
        second.append(mv(ji,[-x for x in forcing],1)+[Adic(0)]*4)
    def curve_dir(v):
        return [[v[0],v[1]],[v[8],v[2]],[v[9],v[3]],[v[10],v[4],v[5],v[6]]]
    columns=[curve_dir(v) for v in directions]
    m=[[col[i] for col in columns] for i in range(4)]
    J=det_poly(m,1);Jdot=[Adic(0)]
    for j in range(4):
        mm=[[v[:] for v in row] for row in m];cc=curve_dir(second[j])
        for i in range(4):mm[i][j]=cc[i]
        Jdot=add(Jdot,det_poly(mm,1),1)
    expected_J=1/(6*xi**5);expected_Jdot=-5/(36*xi**11)
    shortcut_checks=add(J,[-expected_J],1)+add(Jdot,[-expected_Jdot],1)
    assert all(v==0 for v in shortcut_checks)
    Y=[transformed[0],transformed[1]]
    Z=[transformed[10],transformed[4],transformed[5],transformed[6]]
    Zd=curve_dir(eps)[3]
    A=scale(power(Y,3,1),k0**3,1)
    Ad=scale(mul(power(Y,2,1),curve_dir(eps)[0],1),3*k0**3,1)
    AJ=mul(A,J,1);AJd=add(mul(Ad,J,1),mul(A,Jdot,1),1)
    R=add(scale(mul(Z,AJd,1),3,1),scale(mul(AJ,Zd,1),-2,1),1)
    def L(v):return add(scale(mul(Z,diff(v),1),3,1),scale(mul(diff(Z),v,1),-2,1),1)
    rr=R[:];snum=[Adic(0)]
    while len(rr)>5:
        deg=len(rr)-3
        coeff=rr[-1]/(3*Z[3]*(deg-2))
        term=[Adic(0)]*deg+[coeff]
        snum=add(snum,term,1);rr=add(rr,scale(L(term),-1,1),1)
    cols=[L([Adic(0)]*j+[Adic(1)]) for j in range(3)]+[Z,[Adic(0)]+Z]
    mat=[[col[i] if i<len(col) else Adic(0) for col in cols] for i in range(5)]
    sol=mv(inv_adic(mat),rr+[Adic(0)]*(5-len(rr)),1)
    snum=add(snum,sol[:3],1)
    residual=add(R,scale(add(L(snum),mul(Z,sol[3:],1),1),-1,1),1)
    assert all(v==0 for v in residual)
    def record(v):return {'valuation':v.val(),'absolute_precision':v.n,'representative':str(v.q),'zero_to_precision':v==0}
    result={'prime':p,'input_precision':n,'anchor':'X0','target_source':'X0^3*Y0^3','coupling':'5/7',
            'incidence_checked_to':min(v.n for v in vals),'jacobian_inverse_checked_to':min(v.n for row in identity for v in row),
            'minimum_inverse_valuation':min(v.val() for row in ji for v in row),
            'J_degree':len(J)-1,'R_degree':len(R)-1,
            'constant_determinant_shortcut_checked_to':min(Adic(v).n for v in shortcut_checks),
            'compact_remainder':record(sol[3]),'residue_remainder':record(sol[4]),
            'normal_form_reconstruction_to':min(v.n for v in residual),
            'scope':'Precision-tracked adic sample; zeros modulo p^N are not generic vanishing proofs.'}
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--prime',type=int,default=101);ap.add_argument('--precision',type=int,default=60)
    ap.add_argument('--output',type=Path,default=Path(__file__).with_name('source-root-101.json'));a=ap.parse_args()
    main(a.prime,a.precision,a.output)
