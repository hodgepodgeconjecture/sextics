"""Bounded elimination probe; a timeout is not a mathematical obstruction."""
import argparse,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--backend',choices=['sympy','flint'],default='sympy');p.add_argument('--seconds',type=int,default=60);p.add_argument('--worker',action='store_true');args=p.parse_args()
if not args.worker:
 start=time.monotonic()
 try:
  r=subprocess.run([sys.executable,__file__,'--worker','--backend',args.backend],capture_output=True,text=True,timeout=args.seconds)
  result=json.loads(r.stdout) if r.returncode==0 else {'status':'worker_failed','returncode':r.returncode,'error':r.stderr[-2000:]}
 except subprocess.TimeoutExpired:result={'status':'timeout','completed_Groebner_basis':False}
 result.update(backend=args.backend,time_limit_seconds=args.seconds,elapsed_seconds=round(time.monotonic()-start,3),scope='Unsaturated system at K=1,L=2. No complete critical algebra or trace is implied by an incomplete run.')
 (ROOT/f'critical-elimination-{args.backend}.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));sys.exit(0)
import resource
resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
import sympy as s
from source_root_chart import build
x=s.symbols('xi a b c d e f z');la=s.symbols('l0:7');xi,a,b,c,d,e,f,z=x;t=s.Symbol('t')
raw,_,_=build('X1');old=[xi,a,b,c,d,e,f,1,1,2,z]
Q=[s.expand(sum(co*s.prod(v**n for v,n in zip(old,ex)) for co,ex in q)) for q in raw]
D=s.discriminant(z+d*t+e*t*t+f*t**3,t);J=s.Matrix(Q).jacobian(x)
fs=Q+[s.diff(D,v)-sum(la[j]*J[j,i] for j in range(7)) for i,v in enumerate(x)]
if args.backend=='sympy':
 G=s.groebner(fs,*(x+la),modulus=13,order='grevlex',method='f5b')
 result={'status':'completed','prime':13,'completed_Groebner_basis':True,'basis_length':len(G.polys),'zero_dimensional':bool(G.is_zero_dimensional)}
else:
 from flint import fmpz_mpoly_ctx,fmpz_mpoly_vec
 ctx=fmpz_mpoly_ctx.get(tuple(map(str,x+la)),'degrevlex')
 polys=[ctx.from_dict({ex:int(co) for ex,co in s.Poly(q,*(x+la)).terms()}) for q in fs]
 G,success=fmpz_mpoly_vec(polys,ctx).buchberger_naive(limits=(1000,50000,4096))
 result={'status':'completed' if success else 'growth_limit','arithmetic':'integer polynomials','completed_Groebner_basis':bool(success),'basis_length':len(G),'limits':[1000,50000,4096]}
print(json.dumps(result))
