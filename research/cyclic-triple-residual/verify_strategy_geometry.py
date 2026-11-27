"""Exact algebra for the strategy reassessment; not a target certificate.

The bundle calculations are determinant/weight checks. Their geometric use
requires the regularity hypotheses explained in strategy-reassessment.md.
No singularity resolution, coniveau theorem, or rational current is supplied.
"""
import hashlib
import itertools
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
record_bytes = (ROOT / 'inputs/cyclic-triple-residual-instance.rec').read_bytes()
record = record_bytes.decode()
repository_record = ROOT.parents[1] / 'data/cyclic-triple-residual/instance.rec'
if repository_record.exists():
    # The archived copy has one extra terminal newline; all record content
    # must otherwise agree, rather than merely matching a duplicated formula.
    assert repository_record.read_text().rstrip() == record.rstrip()
x = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
X0, X1, X2, Y0, Y1, Y2 = x
F = sp.sympify(' '.join(record.split('F:', 1)[1].split('\nA:', 1)[0].splitlines()),
               locals=dict(zip(map(str, x), x)))
A = sp.sympify(record.split('\nA:', 1)[1].splitlines()[0],
               locals=dict(zip(map(str, x), x)))
assert A == X0**3 * Y0**3
assert record.split('\nq:', 1)[1].splitlines()[0].strip() == '1'
k = sp.Rational(5, 7)
z = sp.Symbol('z')
B = X0**3 + X0*Y1*Y2 + X1*Y0*Y1
C = (sum(v**6 for v in [X0, X1, Y0, Y1, Y2])
     + k*(X0**3*X1*Y0*Y1 + X0**2*Y0*Y1*Y2**2 + X0*X1*Y1**3*Y2))
H = z**2 + k*B*z + C
assert sp.expand(H.subs(z, X2**3) - F) == 0
assert sp.expand(4*H - ((2*z+k*B)**2 - (k**2*B**2-4*C))) == 0
ww = sp.Symbol('ww')
assert sp.expand(H.subs(z,ww**3)
                 -(ww**3-z)*(ww**3+k*B+z)-H) == 0

# Binary-cubic discriminant: verify a generating set for GL(2).
s, t, a, b, c, d, lam, mu = sp.symbols('s t a b c d lam mu')
Z = a*s**3 + b*s**2*t + c*s*t**2 + d*t**3
def disc(binary):
    return sp.discriminant(sp.expand(binary).subs(t, 1), s)
Delta = disc(Z)
assert sp.expand(disc(Z.subs(s, s+lam*t))-Delta) == 0
assert sp.expand(disc(Z.subs({s:t, t:s}, simultaneous=True))-Delta) == 0
assert sp.expand(disc(Z.subs({s:lam*s, t:mu*t}, simultaneous=True))
                 -(lam*mu)**6*Delta) == 0
assert {sum(e) for e, _ in sp.Poly(Delta, a, b, c, d).terms()} == {4}
zero_cubic = {a:0, b:0, c:0, d:0}
assert all(sp.diff(Delta, v).subs(zero_cubic) == 0 for v in [a,b,c,d])

# Tot(Sym^3 U*) over Gr(2,5), equations in Sym^6 U*.
def determinant_exponent(m):
    weights = [(j, m-j) for j in range(m+1)]
    sums = tuple(sum(w[i] for w in weights) for i in range(2))
    assert sums[0] == sums[1]
    return sums[0]
v_det = determinant_exponent(3)
w_det = determinant_exponent(6)
canonical_M = -5-v_det+w_det
canonical_root_model = -5-v_det-1+w_det+6
assert (v_det, w_det, canonical_M, canonical_root_model) == (6,21,10,15)
assert 6+4-7 == 3
# Four evaluations determine a cubic. Each value obeys a monic quadratic;
# this matrix is the linear algebra input in the finiteness argument.
evaluation = sp.Matrix([[r**j for j in range(4)] for r in range(4)])
assert evaluation.det() == 12

# Actual smooth line-space seed in the projective chart X1=1, Y1=t.
xi, aa, bb, cc, dd, ee, ff, k0, k2, k3 = sp.symbols(
    'xi aa bb cc dd ee ff k0 k2 k3')
dependent = [xi, aa, bb, cc, dd, ee, ff]
line_sub = {Y0:xi+aa*t, Y1:t, X0:k0+bb*t, X1:1,
            Y2:k2+cc*t, z:k3+dd*t+ee*t**2+ff*t**3}
seed = {xi:sp.I, aa:0, bb:0, cc:sp.I, dd:0, ee:0, ff:0,
        k0:0, k2:0, k3:0}
incidence = sp.Poly(sp.expand(H.subs(line_sub, simultaneous=True)), t)
equations = sp.Matrix([incidence.nth(j) for j in range(7)])
assert equations.subs(seed) == sp.zeros(7,1)
seed_jacobian = equations.jacobian(dependent).subs(seed)
seed_det = sp.simplify(seed_jacobian.det())
assert seed_det != 0
# Thus M is smooth of dimension 3 here. delta^6=disc(Z) has zero
# differential at this point, so its unnormalized hypersurface model is
# singular. The first jet proves Delta is not identically zero on M.
alpha, beta, gamma = sp.symbols('alpha beta gamma')
endpoint_jet = sp.Matrix([alpha, beta, gamma])
dependent_jet = -seed_jacobian.inv()*equations.jacobian([k0,k2,k3]).subs(seed)*endpoint_jet
Z_jet = sp.expand(gamma+dependent_jet[4]*t+dependent_jet[5]*t**2+dependent_jet[6]*t**3)
assert Z_jet == gamma-alpha*t**3
leading_Delta = sp.factor(sp.discriminant(Z_jet,t))
assert leading_Delta == -27*alpha**2*gamma**2

# The parity quotient over the two ratios Y1/Y0, Y2/Y0.
u, v, y, r = sp.symbols('u v y r')
substituted = sp.Poly(sp.expand(F.subs({Y0:y,Y1:u*y,Y2:v*y}, simultaneous=True)), y)
assert all(e[0] % 2 == 0 for e, _ in substituted.terms())
surface = sp.expand(sum(coeff*r**(e[0]//2) for e,coeff in substituted.terms()))
assert sp.expand(surface.subs(r,y**2)-substituted.as_expr()) == 0
surface_vars = [X0,X1,X2,r]
assert {sum(ei*wi for ei,wi in zip(e,[1,1,1,2]))
        for e,_ in sp.Poly(surface,*surface_vars).terms()} == {6}
assert sp.expand(surface.subs({u:0,v:0})
                 -(X0**6+X1**6+X2**6+r**3+k*X0**3*X2**3)) == 0
# At u=v=0 the gradient forces r=X1=0, then either X0=X2=0
# or a 2x2 system of determinant 4-k^2. Hence quasismoothness.
assert 4-k**2 == sp.Rational(171,49)
canonical_surface = 6-sum([1,1,1,2])
assert canonical_surface == 1
surface_top_forms = {str(q):(weight+1)%3 for q,weight in zip([X0,X1,X2],[0,0,1])}
assert list(surface_top_forms.values()).count(1) == 2

# Full diagonal symmetry and the stabilizer of the target residue.
monomials = [e for e,_ in sp.Poly(F,*x).terms()]
group = [w+(0,) for w in itertools.product(range(6),repeat=5)
         if all(sum(ei*wi for ei,wi in zip(e,w+(0,)))%6 == 0 for e in monomials)]
target_character = lambda w:(sum(w)+3*w[0]+3*w[3])%6
kernel = [w for w in group if target_character(w) == 0]
g = (1,5,5,2,2,0)
assert len(group) == 18 and len(kernel) == 6
assert set(kernel) == {tuple(j*wi%6 for wi in g) for j in range(6)}
involutions = [w for w in group if w != (0,)*6 and all(2*wi%6 == 0 for wi in w)]
assert involutions == [(3,3,3,0,0,0)]
# Fix the invariant ratios X1/X2 and Y0/Y1: the residual P3 coordinates
# are X0,X2,Y1,Y2, of g-weights 1,5,2,0. Canonical degree is 2.
fiber_weights = [1,5,2,0]
fiber_names = ['X0','X2','Y1','Y2']
invariant_two_forms = []
for i in range(4):
    for j in range(i,4):
        if (fiber_weights[i]+fiber_weights[j]+sum(fiber_weights))%6 == 0:
            invariant_two_forms.append(fiber_names[i]+'*'+fiber_names[j])
assert invariant_two_forms == ['X2*X2','Y1*Y1']
# At the zero values of both ratios, all mixed gradient terms have
# coefficient sums <= 6*k < 6. The maximum-coordinate argument proves
# smoothness of this ordinary sextic surface, hence of the generic fiber.
pencil_surface = sp.expand(F.subs({X1:0,Y0:0}))
pencil_vars = [X0,X2,Y1,Y2]
mixed = pencil_surface-sum(q**6 for q in pencil_vars)
gradient_bounds = [sum(abs(co) for _,co in sp.Poly(sp.diff(mixed,q),*pencil_vars).terms())
                   for q in pencil_vars]
assert max(gradient_bounds) == 6*k < 6

report = {
    'input_sha256':hashlib.sha256(record_bytes).hexdigest(),
    'checks':{name:'PASS' for name in [
        'exact_target_and_double_cover', 'binary_cubic_discriminant_covariance',
        'bundle_determinants', 'line_seed_and_discriminant_jet',
        'parity_surface_equation_and_weights', 'diagonal_symmetry_and_characters']},
    'line_space':{'expected_dimension':3,'canonical_Plucker_exponent_if_regular':10,
                  'root_model_dualizing_Plucker_exponent_if_lci':15,
                  'seed_incidence_determinant':str(seed_det),
                  'leading_discriminant':str(leading_Delta),
                  'unnormalized_root_model_at_seed':'singular; tangent dimension 4, local dimension 3'},
    'parity_surface':{'equation':str(surface),'canonical_weight':1,
                      'smooth_specialization_gradient_determinant':str(4-k**2),
                      'holomorphic_two_form_deck_characters_mod_3':surface_top_forms},
    'symmetry':{'group_order':18,'target_kernel_order':6,'kernel_generator':g,
                'projective_involutions':involutions,
                'invariant_two_forms_on_two_pencil_fiber':invariant_two_forms,
                'smooth_pencil_specialization_mixed_gradient_bounds':list(map(str,gradient_bounds))},
    'limitations':[
        'No proof of global smoothness, normality, canonical singularities, or a resolution of the actual line/root spaces.',
        'Canonical exponents do not by themselves prove a birational classification for the actual component.',
        'The character audit rejects only the stated vanishing shortcut, not all correspondences.',
        'No scalar rational primitive, protected Lean certificate, or unrestricted nonexistence proof.']}
(ROOT/'strategy-geometry.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
