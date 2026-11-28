"""Five-variable critical equation over F13; not a completed elimination.
Requires python-flint==0.9.0.
"""
import sys,time,json,itertools
from pathlib import Path
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--export-input');args=parser.parse_args()
ROOT=Path(__file__).resolve().parent
import sympy as s
from flint import nmod_mpoly_ctx
from algebraic_curve_model import model
start=time.monotonic();p=13;M=model();a,b,c,h=M['vars'];xi,K,L,z=M['params'];sv=(a,b,c,xi)
ctx=nmod_mpoly_ctx.get(tuple(map(str,sv)),p,'degrevlex');gens=ctx.gens();zero=ctx.constant(0)
def cv(expr):
 return ctx.from_dict({e:int(v.p)*pow(int(v.q),-1,p)%p for e,v in s.Poly(expr.subs({K:1,L:2,h:1}),*sv).terms()})
D=[cv(q) for q in M['D']];B=[cv(q) for q in M['B']];d0=D[0];k=5*pow(7,-1,p)%p
n2=4*d0*D[2]-D[1]**2;n3=8*d0*d0*D[3]-4*d0*D[1]*D[2]+D[1]**3
P=[64*d0**3*D[4]-4*D[1]*n3-n2*n2,64*d0**4*D[5]-n2*n3,256*d0**5*D[6]-n3*n3]
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def mul(x,y):return (x[0]*y[0]+d0*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def scale(x,k):return (k*x[0],k*x[1])
def power(x,n):
 out=(ctx.constant(1),zero)
 for _ in range(n):out=mul(out,x)
 return out
inv=lambda n:pow(n,-1,p)
Z=[(ctx.constant(-k*inv(2)%p),ctx.constant(inv(2))),(-k*inv(2)*B[1]*d0,inv(4)*D[1]),(-k*inv(2)*B[2]*d0*d0,inv(16)*n2),(-k*inv(2)*B[3]*d0**3,inv(32)*n3)]
terms=[(1,[0,2,2,0]),(-4,[0,3,0,1]),(-4,[1,0,3,0]),(-27,[2,0,0,2]),(18,[1,1,1,1])]
U,V=zero,zero
for co,ee in terms:
 q=(ctx.constant(co%p),zero)
 for zz,e in zip(Z,ee):q=mul(q,power(zz,e))
 U,V=add((U,V),q)
q=6
while q and not divmod(U,d0)[1] and not divmod(V,d0)[1]:U=U//d0;V=V//d0;q-=1
print('source numerator built',len(U),len(V),'denominator exponent',q,flush=True)
J=[[f.derivative(j) for j in range(4)] for f in P]
def det(A):
 ans=zero
 for perm in itertools.permutations(range(len(A))):
  term=ctx.constant(((-1)**sum(perm[i]>perm[j] for i in range(len(A)) for j in range(i+1,len(A))))%p)
  for i,j in enumerate(perm):term*=A[i][j]
  ans+=term
 return ans
N=[(-1)**(3+j)*det([[row[i] for i in range(4) if i!=j] for row in J]) for j in range(4)]
print('minors built',[len(f) for f in N],flush=True)
C0,C1=zero,zero
for i in range(4):
 dp=d0.derivative(i)
 C0+=N[i]*(d0*U.derivative(i)-q*dp*U)
 C1+=N[i]*(d0*V.derivative(i)+(inv(2)-q)*dp*V)
print('critical equation built',len(C0),len(C1),flush=True)
removed=0
while not divmod(C0,d0)[1] and not divmod(C1,d0)[1]:
 C0=C0//d0;C1=C1//d0;removed+=1
print('removed D0 power',removed,'terms',len(C0),len(C1),flush=True)
before=[len(C0),len(C1)]
for cycle in range(3):
 oldpair=(C0,C1)
 new=[]
 for cc in oldpair:
  for ff in P:
   quotient,remainder=divmod(cc,ff)
   assert cc==quotient*ff+remainder
   cc=remainder
  new.append(cc)
 C0,C1=new
 print('reduction cycle',cycle,'terms',len(C0),len(C1),flush=True)
 if (C0,C1)==oldpair:break
def ev(f,pt):return sum(int(co)*s.prod(pow(v,e,p) for v,e in zip(pt,ex)) for ex,co in f.to_dict().items())%p
checks=[]
for pt,w,expected in [([6,2,12,4],12,0),([10,12,3,9],12,4)]:
 assert all(ev(f,pt)==0 for f in P)
 numerator=(ev(C0,pt)+w*ev(C1,pt))%p
 denom=pow(int(ev(d0,pt)),q+1-removed,p)*int(ev(N[3],pt))%p
 assert denom and numerator*pow(denom,-1,p)%p==expected
 checks.append({'point_a_b_c_xi':pt,'w':w,'Delta_xi':expected})
report={'prime':p,'K':1,'L':2,'slope_equation_term_counts':[len(f) for f in P],'slope_equation_degrees':[f.total_degree() for f in P], 'Delta_numerator_terms':[len(U),len(V)],'Delta_denominator_D0_exponent':q,'removed_critical_D0_power':removed,'critical_terms_before_incidence_reduction':before,'critical_terms':[len(C0),len(C1)],'critical_degrees':[C0.total_degree(),C1.total_degree()],'checks':checks,'elapsed_seconds':time.monotonic()-start}
(ROOT/'reduced-critical.json').write_text(json.dumps(report,indent=2,default=int)+'\n')
# Solver input: five equations, without saturation. Preserve exact sparse expressions.
if args.export_input:
 ctx5=nmod_mpoly_ctx.get(('a','b','c','xi','w'),p,'degrevlex')
 def embed(f):return ctx5.from_dict({e+(0,):int(co) for e,co in f.to_dict().items()})
 critical5=ctx5.from_dict({**{e+(0,):int(co) for e,co in C0.to_dict().items()},**{e+(1,):int(co) for e,co in C1.to_dict().items()}})
 exported=[embed(d0)-ctx5.gens()[4]**2]+[embed(f) for f in P]+[critical5]
 Path(args.export_input).write_text('a,b,c,xi,w\n13\n'+',\n'.join(map(str,exported))+'\n')
print(json.dumps(report,indent=2,default=int))
