"""Exact modular factor constraints for the generic incidence field.

The characteristic-zero specialization and smooth-seed arguments are in
incidence-irreducibility.md. This does not prove Kummer irreducibility or
rational exactness of the source.
"""
import json
from pathlib import Path
import sympy as s
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_pow_mod,gf_gcd,gf_sub,gf_mul
from finite_curve_algebra import run
from algebraic_curve_model import model,finite_test

ROOT=Path(__file__).resolve().parent


def rabin(f,p):
    """Rabin's irreducibility criterion, with monic descending coefficients."""
    n=len(f)-1
    assert n>0 and f[0]==1
    X=[1,0]
    if n==1:return
    assert not gf_sub(gf_pow_mod(X,p**n,f,p,ZZ),X,p,ZZ)
    for q in s.factorint(n):
        g=gf_sub(gf_pow_mod(X,p**(n//int(q)),f,p,ZZ),X,p,ZZ)
        assert gf_gcd(f,g,p,ZZ)==[1]


def seed_check():
    M=model();a,b,c,h=M['vars'];xi,K,L,z0=M['params']
    t,d,e,f=s.symbols('t d e f');kap=s.Rational(5,7)
    B=sum(q*t**j for j,q in enumerate(M['B']))
    D=sum(q*t**j for j,q in enumerate(M['D']))
    Z=z0+d*t+e*t*t+f*t**3
    H=s.Poly(s.expand(Z*Z+kap*B*Z+(kap*kap*B*B-D)/4).subs(h,1),t)
    Q=s.Matrix([H.nth(j) for j in range(7)])
    dep=[xi,a,b,c,d,e,f]
    seed={xi:s.I,a:0,b:0,c:s.I,d:0,e:0,f:0,K:0,L:0,z0:0}
    assert Q.subs(seed)==s.zeros(7,1)
    J=Q.jacobian(dep).subs(seed);det=s.simplify(J.det())
    assert det==-135000*s.I/2401
    aa,bb,cc=s.symbols('alpha beta gamma')
    jet=-J.inv()*Q.jacobian([K,L,z0]).subs(seed)*s.Matrix([aa,bb,cc])
    zjet=s.expand(cc+jet[4]*t+jet[5]*t*t+jet[6]*t**3)
    assert zjet==cc-aa*t**3
    assert s.factor(s.discriminant(zjet,t))==-27*aa*aa*cc*cc
    return {'field':'Q(i)','incidence_Jacobian_determinant':str(det),
            'endpoint_w_derivative_in_k3':2,'first_cubic_jet':str(zjet),
            'smooth_rational_seed_and_dominance':'PASS'}


def main():
    X=s.Symbol('T');records=[];possible=set(range(1,120))
    for p in [1009,1013,1069]:
        assert s.isprime(p) and p%4==1
        fiber=run(p)
        infinity=finite_test(model(),p)
        assert infinity['endpoint']==fiber['endpoint'] and infinity['rank']==105
        f=fiber['discriminant_minimal_polynomial_coefficients_ascending'][::-1]
        rec={key:fiber[key] for key in ['prime','endpoint','w0_mod_p','macaulay_ranks','quotient_dimension','checks']}
        rec['minimal_polynomial_coefficients_ascending']=f[::-1]
        rec['no_direction_infinity_certificate']={
            'macaulay_degree':13,'matrix_shape':infinity['matrix_shape'],
            'rank':105,'selected_independent_rows':infinity['selected_independent_rows']}
        if p!=1069:
            poly=s.Poly.from_list(f,X,modulus=p)
            fac=s.factor_list(poly)[1];factors=[];product=[1]
            for q,multiplicity in fac:
                assert multiplicity==1
                coeff=[int(v)%p for v in q.all_coeffs()]
                rabin(coeff,p)
                product=gf_mul(product,coeff,p,ZZ)
                factors.append(coeff)
            assert product==f
            degrees=sorted(len(q)-1 for q in factors)
            sums={0}
            for n in degrees:sums|={v+n for v in list(sums)}
            possible&=sums
            rec.update({'irreducible_factor_degrees':degrees,
                        'irreducible_factors_coefficients_descending':factors,
                        'rabin_and_product_checks':'PASS',
                        'possible_proper_factor_degrees_after_this_fiber':sorted(possible)})
        else:
            rem=gf_sub(gf_pow_mod([1,0],p,f,p,ZZ),[1,0],p,ZZ)
            assert gf_gcd(f,rem,p,ZZ)==[1]
            assert possible=={1,119}
            possible-={1,119}
            rec.update({'gcd_with_T_to_p_minus_T':1,'linear_factor_exclusion':'PASS',
                        'possible_proper_factor_degrees_after_this_fiber':sorted(possible)})
        records.append(rec)
    assert not possible
    report={'fibers':records,'seed':seed_check(),
            'conclusion':'Generic degree-120 incidence polynomial irreducible over Q(i)(K,L,xi,w0); smooth rational seed yields a geometrically integral selected incidence component. See note for specialization argument.',
            'limitations':['Not irreducibility of every specialized endpoint curve.',
                           'Not irreducibility of the sixth-root discriminant extension.',
                           'No rational primitive, target certificate, or Lean formalization.']}
    # Together with the degree-two boundary point, this also excludes an
    # algebraic constant extension in the generic curve over C(K,L).
    p=records[0]['prime'];value=-367*pow(49,-1,p)%p
    assert 1 in records[0]['irreducible_factor_degrees']
    assert pow(value,(p-1)//2,p)==p-1
    report['relative_constant_field_test']={
        'candidate_extension':'sqrt((kappa^2-4)*K^6-4)',
        'candidate_square_at_seed':-4,'seed_square_roots_in_Q_i':['2*i','-2*i'],
        'prime':p,'value_at_K_one':value,'Euler_criterion':p-1,
        'rational_incidence_point_exists':True,
        'scope':'Use the smooth degree-two boundary point from finite-boundary.md and the descent argument in incidence-irreducibility.md.'}
    (ROOT/'incidence-irreducibility.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps([{'prime':r['prime'],'factors':r.get('irreducible_factor_degrees'),
                      'remaining_degrees':r['possible_proper_factor_degrees_after_this_fiber']} for r in records],indent=2))


if __name__=='__main__':main()
