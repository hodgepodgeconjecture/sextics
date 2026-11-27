"""A transverse, horizontal discriminant zero of the fixed incidence family.

The small finite-field witness lifts in characteristic zero by the Jacobian
criterion. It certifies genuine degree-six Kummer branching, not exactness.
"""
import json
from pathlib import Path
import sympy as s
from algebraic_curve_model import model
from line_exactness import inverse,mv

ROOT=Path(__file__).resolve().parent


def run():
    M=model();a,b,c,h=M['vars'];xi,K,L,z0=M['params']
    t,d,e,f=s.symbols('t d e f');kap=s.Rational(5,7)
    B=sum(q*t**j for j,q in enumerate(M['B']))
    D=sum(q*t**j for j,q in enumerate(M['D']))
    Z=z0+d*t+e*t*t+f*t**3
    H=s.Poly(s.expand(7*(Z*Z+kap*B*Z+(kap*kap*B*B-D)/4)).subs(h,1),t)
    Q=s.Matrix([H.nth(j) for j in range(7)])
    variables=[xi,a,b,c,d,e,f,K,L,z0]
    values=[9,10,12,3,6,8,3,1,2,1];p=13
    point=dict(zip(variables,values))
    assert all(int(q.subs(point))%p==0 for q in Q)
    J=Q.jacobian(variables)
    full=[[int(q.subs(point))%p for q in row] for row in J.tolist()]
    dep=[row[:7] for row in full];ji=inverse(dep,p)
    det=int(s.Matrix(dep).det())%p;assert det
    Delta=s.discriminant(Z,t)
    assert int(Delta.subs(point))%p==0
    dg=[int(Delta.diff(v).subs(point))%p for v in variables]
    derivatives=[]
    for j in [7,8,9]:
        dy=mv(ji,[-row[j] for row in full],p)
        velocity=dy+[int(q==j) for q in [7,8,9]]
        derivatives.append(sum(x*y for x,y in zip(dg,velocity))%p)
    assert derivatives==[9,3,6]
    augmented=[row[:7]+[row[9]] for row in full]+[dg[:7]+[dg[9]]]
    augdet=int(s.Matrix(augmented).det())%p;assert augdet==det*derivatives[2]%p
    cubic=s.Poly(Z.subs(point),t,modulus=p)
    assert cubic.degree()==3 and s.gcd(cubic,cubic.diff()).degree()==1
    # True source deformation, keeping xi,K,L fixed and allowing z0 to move.
    A=s.Poly(sum(q*t**j for j,q in enumerate(M['A'])).subs(h,1),t)
    active=[1,2,3,4,5,6,9]
    js=[[row[j] for j in active] for row in full]
    rhs=[int((7*A.nth(j)).subs(point))%p for j in range(7)]
    ds=mv(inverse(js,p),rhs,p)
    dot=sum(dg[j]*q for j,q in zip(active,ds))%p
    k=5*pow(7,-1,p)%p;w=(2*point[z0]+k*point[K]**3)%p
    xi_k3=-w*pow(6*pow(point[xi],5,p),-1,p)%p
    delta_xi=derivatives[2]*pow(xi_k3,-1,p)%p
    assert ds[-1]==point[K]**3*point[xi]**3*pow(w,-1,p)%p
    pairing=-2*dot*pow(w*w*delta_xi*delta_xi,-1,p)%p
    return {'prime':p,'variables':list(map(str,variables)),'point':values,
            'incidence_Jacobian_determinant_mod_p':det,
            'augmented_incidence_discriminant_Jacobian_determinant_mod_p':augdet,
            'discriminant_endpoint_derivatives_k0_k2_k3':derivatives,
            'cubic':str(cubic.as_expr()),'gcd_with_cubic_derivative_degree':1,
            'w0':w,'xi_k3':xi_k3,'Delta_xi':delta_xi,
            'Delta_source_derivative_at_fixed_xi':dot,
            'local_opposite_character_pairing_residue':pairing,
            'checks':{'fixed_target_binding':'PASS','incidence_point':'PASS',
                      'smooth_incidence':'PASS','horizontal_simple_discriminant_zero':'PASS',
                      'ordinary_cubic_double_root':'PASS','actual_source_variation':'PASS'},
            'scope':'The horizontal simple zero proves generic degree-six Kummer branching. The local pairing contribution is not the global obstruction sum. No target certificate.'}


if __name__=='__main__':
    report=run()
    (ROOT/'kummer-branch.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
