"""Exact local identities for a conditional global principal-parts test.
No global critical-locus trace or rational primitive is computed here.
"""
import json
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parent
x,K,w,D,Dx,Ds,Dxx,Dsx,r=s.symbols('x K w D Dx Ds Dxx Dsx r',nonzero=True)
wprime=-12*x**5/w
V=2*K**3*x**3/(3*w**3)+Ds/(18*w*D)
g=-Ds/(3*w*Dx)
gprime=s.diff(g,Ds)*Dsx+s.diff(g,w)*wprime+s.diff(g,Dx)*Dxx
R=s.cancel(V-(gprime-g*Dx/(6*D)))
assert s.cancel(R-(2*K**3*x**3/(3*w**3)-gprime))==0
# Simple discriminant zero: Delta=z^6, dxi=6*z^5/Dx dz.
z=s.symbols('z')
alpha=Ds/(18*w*z**7)*6*z**5/Dx
primitive=-Ds/(3*w*Dx*z)
assert s.diff(primitive,z)==alpha
theta=6/(w*Dx)
simple=s.residue(primitive*theta,z,0)
assert s.cancel(simple+2*Ds/(w*w*Dx*Dx))==0
# Simple critical point: Delta_xi=Dxx*u+O(u^2).
u=s.symbols('u')
critical=6*s.residue(-Ds/(3*w*w*D*Dxx*u),u,0)
assert critical==-2*Ds/(w*w*D*Dxx)
# Endpoint w=0: dxi=-w/(12*xi^5) dw.
endpoint=6*s.residue(K**3/(18*x*x*u*u*D)*(-u/(12*x**5)),u,0)
assert endpoint==-K**3/(36*x**7*D)
T,Ns,LN=s.symbols("T Ns LN",nonzero=True)
critical_second=36*x**10*LN/(w*w*T*T)
critical_trace=2*(Ns/T)/(w*w*D*critical_second)
assert s.cancel(critical_trace-T*Ns/(18*x**10*D*LN))==0
orders=[]
for name,e,m,ell,at_infinity in [
 ('infinity excess family',1,-12,0,True),
 ('infinity c=0',4,-24,3,True),
 ('infinity b=0',4,-20,7,True),
 ('finite constant cubic',5,-18,3,False),
 ('finite translated cubic',4,-24,4,False)]:
 order=(4*e-ell-m-1) if at_infinity else (2*e-ell-m-1)
 assert order>=0
 orders.append({'branch':name,'chi_order_lower_bound':order})
# The previous actual-source witness checks the simple-zero formula.
p=13
assert (-2*4*pow(12**2*4**2,-1,p))%p==6
report={'checks':{'exact_discriminant_pole_subtraction':'PASS','simple_zero_pairing':'PASS',
 'critical_point_pairing':'PASS','endpoint_pairing':'PASS','audited_boundary_chi_regularity':'PASS',
 'previous_actual_source_witness':'PASS'},
 'regularized_density':str(R),'simple_zero_pairing':str(simple),
 'six_sheet_critical_residue':str(critical),'six_sheet_endpoint_residue':str(endpoint),
 'boundary_orders':orders,'cofactor_critical_trace_summand':str(s.factor(critical_trace)),
 'scope':'Conditional residue-duality reduction. Remaining boundary terms and global regularity of theta are not established; the critical trace has not been evaluated.'}
(ROOT/'principal-parts.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
