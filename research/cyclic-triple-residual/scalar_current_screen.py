"""Precision-tracked screen of endpoint-polynomial scalar currents.

Ansatz: b_i=J*M_i/xi^m, where M_i is a constant-coefficient weighted
homogeneous polynomial of degree weight(k_i)+m in (xi,k0,k1,k2,k3).
This is a bounded class, not all rational currents in the incidence field.
Euler currents remove a known kernel before testing the remaining columns.
"""
import argparse
import itertools
import hashlib
import json
import random
import time
from pathlib import Path
from fractions import Fraction
import sympy as sp
from line_exactness import (compile_poly, derivative, evaluate, inverse, mv,
                            mul, power)
from source_root_chart import Adic, build

WEIGHTS = [1,1,1,1,3]
ENDPOINT_WEIGHTS = [1,1,1,3]


def verify_cleared_formula():
    x,M,Md,xd,Delta,Delta_d,Delta_s,k0=sp.symbols('x M Md xd Delta Delta_d Delta_s k0',nonzero=True)
    J=1/(6*x**5)
    Js=-k0**3/(18*x**8)
    for m in [0,2]:
        b=J*M/x**m
        derivative=sp.diff(b,M)*Md+sp.diff(b,x)*xd
        cleared=18*Delta*x**(m+1)*(derivative-b*Delta_d/(6*Delta))/J
        expected=18*Delta*x*Md-18*(m+5)*Delta*M*xd-3*Delta_d*x*M
        assert sp.cancel(cleared-expected)==0
        rhs=18*Delta*x**(m+1)*(Js/3-J*Delta_s/(18*Delta))/J
        assert sp.cancel(rhs-(-2*Delta*k0**3*x**(m-2)-x**(m+1)*Delta_s))==0


def monomials(degree):
    return [e for e in itertools.product(range(degree+1), repeat=5)
            if sum(a*b for a,b in zip(e, WEIGHTS)) == degree]


def columns(m, linear_slopes=False, mixed=False):
    # The exact zero current f*S, f=M/xi^m of weight zero, permits setting
    # every k0 component divisible by k0 to zero. Multiplication by k0
    # gives distinct monomials, so this removes exactly one full Euler copy.
    slopes = [-1,1,2,3,6] if linear_slopes else [-1]
    all_cols = [(i,e,l) for i,w in enumerate(ENDPOINT_WEIGHTS)
                for e in monomials(w+m) for l in slopes]
    if mixed:
        assert m == 2 and not linear_slopes
        for i,w in enumerate(ENDPOINT_WEIGHTS):
            for e in monomials(w):
                lifted=(e[0]+2,*e[1:])
                all_cols.extend((i,lifted,l) for l in [1,2,3,6])
    removed = [(i,e,l) for i,e,l in all_cols if i == 0 and e[1] > 0]
    assert len(removed) == len(monomials(m))*len(slopes)+(4 if mixed else 0)
    return [c for c in all_cols if c not in removed], len(all_cols), len(removed)


def scalar_row(p, n, endpoint_direction, data, cols, m):
    Q, D, disc_terms, disc_derivatives = data
    modulus = p**n
    root_i = next(i for i in range(p) if (i*i+1) % p == 0)
    alpha,beta,gamma = endpoint_direction
    point = [root_i,0,0,root_i,0,0,0,p*alpha,1,p*beta,p*gamma]
    for _ in range(n.bit_length()+2):
        jac = [[evaluate(v,point,modulus) for v in row[:7]] for row in D]
        correction = mv(inverse(jac,modulus), [evaluate(q,point,modulus) for q in Q], modulus)
        point[:7] = [(x-v) % modulus for x,v in zip(point[:7],correction)]
    assert all(evaluate(q,point,modulus) == 0 for q in Q)
    jac = [[evaluate(v,point,modulus) for v in row] for row in D]
    ji = inverse([row[:7] for row in jac],modulus)
    dirs = [mv(ji,[-row[j] for row in jac],modulus)
            +[int(j == l) for l in range(7,11)] for j in range(7,11)]
    xi,a,b,c,d,e,f,k0,k1,k2,k3 = point
    Z = [f,e,d,k3]
    Delta = evaluate(disc_terms,Z,modulus)
    grad = [evaluate(v,Z,modulus) for v in disc_derivatives]
    Delta_i = [sum(g*v for g,v in zip(grad,[u[6],u[5],u[4],u[10]])) % modulus for u in dirs]
    assert (sum(w*k*v for w,k,v in zip(ENDPOINT_WEIGHTS,point[7:],Delta_i))-6*Delta) % modulus == 0
    assert (sum(w*k*u[0] for w,k,u in zip(ENDPOINT_WEIGHTS,point[7:],dirs))-xi) % modulus == 0
    for j in [1,2,3,6]:
        assert sum(w*k*u[j] for w,k,u in zip(ENDPOINT_WEIGHTS,point[7:],dirs)) % modulus == 0
    A = mul(power([k0,b],3,modulus),power([xi,a],3,modulus),modulus)
    A += [0]*(7-len(A))
    U = mv(ji,[7*v for v in A],modulus)
    Delta_s = sum(g*v for g,v in zip(grad,[U[6],U[5],U[4],0])) % modulus
    endpoint = [xi,k0,k1,k2,k3]
    def monomial(exps):
        out = 1
        for x,e in zip(endpoint,exps): out = out*pow(x,e,modulus) % modulus
        return out
    def partial(exps,j):
        if not exps[j]: return 0
        e = list(exps); e[j] -= 1
        return exps[j]*monomial(e) % modulus
    values = []
    for i,exps,l in cols:
        base = monomial(exps)
        base_d = (partial(exps,0)*dirs[i][0]+partial(exps,i+1)) % modulus
        M = base if l == -1 else base*point[l] % modulus
        Md = base_d if l == -1 else (base_d*point[l]+base*dirs[i][l]) % modulus
        # 18*Delta*xi^(m+1) times the equation divided by J.
        col = (18*Delta*xi*Md - 18*(m+5)*Delta*M*dirs[i][0]
               -3*Delta_i[i]*xi*M) % modulus
        values.append(col)
    # B/J=-k0^3/(9*xi^3)-Delta_s/(18*Delta).
    rhs = (-2*Delta*pow(k0,3,modulus)*pow(xi,m-2,modulus)
           -pow(xi,m+1,modulus)*Delta_s) % modulus
    return [Adic(v,n) for v in values+[rhs]]


def screen(p,n,m,sample_limit,seed,linear_slopes=False,mixed=False):
    start = time.monotonic()
    verify_cleared_formula()
    Adic.p = p; Adic.cap = 2*n
    Q,D,_ = build('X1')
    t,a,b,c,d = sp.symbols('t a b c d')
    dc = compile_poly(sp.discriminant(a*t**3+b*t*t+c*t+d,t),[a,b,c,d])
    data = (Q,D,dc,[derivative(dc,j) for j in range(4)])
    cols, original_count, euler_count = columns(m,linear_slopes,mixed)
    rng = random.Random(seed)
    basis = []
    used = []
    unpivoted = set(range(len(cols)))
    witness = None
    pivot_records = []
    for sample in range(sample_limit):
        direction = [rng.randrange(1,p) for _ in range(3)]
        row = scalar_row(p,n,direction,data,cols,m)
        used.append(direction)
        for pivot,old in basis:
            factor = row[pivot]
            # A zero representative still carries an uncertainty O(p^N).
            # Propagate it instead of silently treating it as an exact zero.
            row = [x-factor*y for x,y in zip(row,old)]
            row[pivot] = Adic(0)  # Exact zero from the row operation.
        resolved = [j for j in unpivoted if row[j] != 0]
        if resolved:
            pivot = min(resolved,key=lambda j:row[j].val())
            # An unresolved coefficient cannot be treated as exact zero.
            assert min(row[j].n for j in unpivoted) > row[pivot].val(), 'increase precision'
            val = row[pivot].val()
            inv = row[pivot].inv()
            row = [v*inv for v in row]
            row[pivot] = Adic(1)
            basis.append((pivot,row))
            unpivoted.remove(pivot)
            pivot_records.append({'sample':sample,'column':pivot,'valuation':val,
                                  'minimum_precision':min(v.n for v in row)})
        elif not unpivoted and row[-1] != 0:
            witness = {'sample':sample, 'rhs_valuation':row[-1].val(),
                       'rhs_precision':row[-1].n,
                       'normalized_rhs_mod_p':int(row[-1].q/(Fraction(p)**row[-1].val())) % p}
            break
        if (sample+1) % 20 == 0:
            print(f'samples={sample+1} rank={len(basis)}/{len(cols)} elapsed={time.monotonic()-start:.1f}s',flush=True)
    result = {'prime':p,'input_precision':n,'denominator_exponent':m,'seed':seed,
              'coefficient_span':['1','a','b','c','f'] if linear_slopes else ['1'],
              'additional_polynomial_slope_span':['a','b','c','f'] if mixed else [],
              'original_columns':original_count,'exact_Euler_kernel_dimension':euler_count,
              'columns_after_Euler_gauge':len(cols),'rank_lower_bound':len(basis),
              'sample_directions':used,'pivots':pivot_records,'inconsistency_witness':witness,
              'status':'EXCLUDED' if witness else 'INCONCLUSIVE',
              'scope':'Only the stated finite polynomial coefficient span for b_i=J*M_i/xi^m. Does not exclude other incidence-dependent currents or general rational denominators.',
              'elapsed_seconds':round(time.monotonic()-start,3)}
    return result


def compact_report(result):
    result=dict(result)
    directions=result.pop('sample_directions')
    pivots=result.pop('pivots')
    def digest(data):
        return hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    result['samples_checked']=len(directions)
    result['reproduction_hashes']={'sample_directions_sha256':digest(directions),
                                  'pivot_records_sha256':digest(pivots)}
    result['pivot_statistics']={'minimum_precision':min(p['minimum_precision'] for p in pivots),
                                'maximum_valuation':max(p['valuation'] for p in pivots),
                                'sum_of_valuations':sum(p['valuation'] for p in pivots)}
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--prime',type=int,default=101)
    ap.add_argument('--precision',type=int,default=36)
    ap.add_argument('--denominator-exponent',type=int,default=2)
    ap.add_argument('--sample-limit',type=int,default=145)
    ap.add_argument('--seed',type=int,default=20261126)
    group=ap.add_mutually_exclusive_group()
    group.add_argument('--linear-slopes',action='store_true')
    group.add_argument('--mixed',action='store_true')
    ap.add_argument('--full-report',action='store_true',help='retain all reproducible directions and pivot records')
    ap.add_argument('--output',type=Path,default=Path(__file__).with_name('scalar-screen-101.json'))
    args=ap.parse_args()
    result=screen(args.prime,args.precision,args.denominator_exponent,args.sample_limit,args.seed,args.linear_slopes,args.mixed)
    if not args.full_report: result=compact_report(result)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['sample_directions','pivots']},indent=2))
