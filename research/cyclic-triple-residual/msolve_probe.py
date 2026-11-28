"""Export and run the critical multiplier system with a bounded msolve call."""
import argparse,json,os,resource,subprocess,tempfile,time
from pathlib import Path
import sympy as s
from source_root_chart import build
from algebraic_curve_model import model
p=argparse.ArgumentParser();p.add_argument('--msolve',required=True);p.add_argument('--prime',type=int,default=1000003);p.add_argument('--seconds',type=int,default=60);p.add_argument('--threads',type=int,default=2);p.add_argument('--memory-gib',type=int,default=3);p.add_argument('--saturate',action='store_true');p.add_argument('--output',required=True);args=p.parse_args()
assert s.isprime(args.prime) and args.prime not in [2,3,7]
model()
x=s.symbols('xi a b c d e f z');la=s.symbols('l0:7');xi,a,b,c,d,e,f,z=x;t=s.Symbol('t')
raw,_,_=build('X1');old=[xi,a,b,c,d,e,f,1,1,2,z]
Q=[s.expand(sum(co*s.prod(v**n for v,n in zip(old,ex)) for co,ex in q)) for q in raw]
D=s.discriminant(z+d*t+e*t*t+f*t**3,t);J=s.Matrix(Q).jacobian(x)
fs=Q+[s.diff(D,v)-sum(la[j]*J[j,i] for j in range(7)) for i,v in enumerate(x)]
variables=x+la
if args.saturate:
 inv=s.Symbol('inv');variables=(inv,)+variables;fs.append(s.expand(inv*xi*(14*z+5)*D-1))
out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
inp=out/'input.ms';inp.write_text(','.join(map(str,variables))+'\n'+str(args.prime)+'\n'+',\n'.join(str(s.Poly(q,*variables).as_expr()).replace('**','^') for q in fs)+'\n')
cmd=[args.msolve,'-g','1','-l','2','-t',str(args.threads),'-v','2','-m','1000','-f',str(inp),'-o',str(out/'basis.out')]
# Use the ordinary F4 engine with a Rabinowitsch variable, not native -S.
def limit():resource.setrlimit(resource.RLIMIT_AS,(args.memory_gib*1024**3,args.memory_gib*1024**3))
env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';start=time.monotonic()
with (out/'solver.log').open('w') as log:
 try:
  result=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=args.seconds,preexec_fn=limit)
  status='completed' if result.returncode==0 else 'solver_error';code=result.returncode
 except subprocess.TimeoutExpired:status='timeout';code=None
report={'prime':args.prime,'parameters':[1,2],'saturation':'xi*(14*z+5)*Delta' if args.saturate else None,'threads':args.threads,'timeout_seconds':args.seconds,'address_space_GiB':args.memory_gib,'status':status,'exit_code':code,'elapsed_seconds':time.monotonic()-start,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'scope':'No trace computed by this runner; saturation does not yet exclude every exceptional stratum.'}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
