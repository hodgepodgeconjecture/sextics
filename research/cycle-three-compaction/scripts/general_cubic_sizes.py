"""Covariant sizes for the general ten-coefficient ternary cubic, using the
SageMath TernaryCubic definitions (Hessian, Theta via polar-conic adjugates,
J = Jacobian determinant / 9). Restricting to the cyclic family reproduces
the 14 / 207 / 1095 terms of data/cycle-three/CubicCovariants.lean."""
import flint, time
names=('a30','a03','a00','a21','a20','a12','a02','a10','a01','a11','x','y','z')
ctx=flint.fmpq_mpoly_ctx.get(names,'lex'); g=ctx.gens()
a30,a03,a00,a21,a20,a12,a02,a10,a01,a11,x,y,z=g
U=a30*x**3+a03*y**3+a00*z**3+a21*x**2*y+a20*x**2*z+a12*x*y**2+a02*y**2*z+a10*x*z**2+a01*y*z**2+a11*x*y*z
def D(p,v): return p.derivative(names.index(v))
def det3(M):
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
V=['x','y','z']
H=det3([[D(D(U,i),j) for j in V] for i in V])*flint.fmpq(1,216)
def polar_adj(P):
    # coefficients of P as cubic in x,y,z, Sage polar_conic then adjugate
    def co(m):
        return P.coeff(m) if False else None
    return None
# Build Theta via Sage recipe: polar conic matrix A(x) with entries from coefficients of the cubic form
def coeffs_of(P):
    # P homogeneous cubic in x,y,z with coefficients polynomials in a's
    mons={'x3':(3,0,0),'y3':(0,3,0),'z3':(0,0,3),'x2y':(2,1,0),'x2z':(2,0,1),'xy2':(1,2,0),'y2z':(0,2,1),'xz2':(1,0,2),'yz2':(0,1,2),'xyz':(1,1,1)}
    out={k:ctx.from_dict({}) for k in mons}
    for exps,cf in P.to_dict().items():
        e=exps[10:]
        for k,m in mons.items():
            if tuple(e)==m:
                out[k]+=ctx.from_dict({tuple(list(exps[:10])+[0,0,0]):cf})
    return out
def polar_matrix(P):
    c=coeffs_of(P); h=flint.fmpq(1,2)
    A00=3*x*c['x3']+y*c['x2y']+z*c['x2z']; A11=x*c['xy2']+3*y*c['y3']+z*c['y2z']; A22=x*c['xz2']+y*c['yz2']+3*z*c['z3']
    A01=x*c['x2y']+y*c['xy2']+h*z*c['xyz']; A02=x*c['x2z']+h*y*c['xyz']+z*c['xz2']; A12=h*x*c['xyz']+y*c['y2z']+z*c['yz2']
    return [[A00,A01,A02],[A01,A11,A12],[A02,A12,A22]]
def adj(M):
    C=[[None]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            r=[k for k in range(3) if k!=i]; s=[k for k in range(3) if k!=j]
            C[j][i]=(-1)**(i+j)*(M[r[0]][s[0]]*M[r[1]][s[1]]-M[r[0]][s[1]]*M[r[1]][s[0]])
    return C
def conic_coeffs(A): return (A[0][0],A[1][1],A[2][2],A[0][1],A[0][2],A[1][2])
t=time.time()
Ua=conic_coeffs(adj(polar_matrix(U))); Ha=conic_coeffs(adj(polar_matrix(H)))
a0,b0,c0,h0,g0,f0=Ua; a1,b1,c1,h1,g1,f1=Ha
Th=flint.fmpq(1,9)*((b0*c1+c0*b1-2*f0*f1)*x**2+(a0*c1+c0*a1-2*g0*g1)*y**2+(a0*b1+b0*a1-2*h0*h1)*z**2+2*(f0*g1+g0*f1-c0*h1-h0*c1)*x*y+2*(h0*f1+f0*h1-b0*g1-g0*b1)*x*z+2*(g0*h1+h0*g1-a0*f1-f0*a1)*y*z)
J=flint.fmpq(1,9)*det3([[D(P,v) for v in V] for P in (U,H,Th)])
print('general cubic: U',len(U),'H',len(H),'Theta',len(Th),'J',len(J), f'{time.time()-t:.1f}s')
# restrict to cyclic family to compare normalization with PR (a21=d,a02=e,a10=f, others 0 besides a30,a03,a00)
sub={}
import itertools
def restrict(P):
    d={}
    for exps,cf in P.to_dict().items():
        # zero out a20,a12,a01,a11 (indices 4,5,8,9)
        if exps[4] or exps[5] or exps[8] or exps[9]: continue
        d[exps]=cf
    return ctx.from_dict(d)
print('restricted to cyclic family: H',len(restrict(H)),'Theta',len(restrict(Th)),'J',len(restrict(J)))
