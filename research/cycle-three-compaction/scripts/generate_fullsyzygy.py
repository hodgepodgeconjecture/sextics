"""Generate FullSyzygy.lean (the whole cyclic covariant syzygy as one Kronecker
integer identity) from the covariant formulas of data/cycle-three/CubicCovariants.lean.

usage: python3 generate_fullsyzygy.py <repo-root> > FullSyzygy.lean
"""
import re, sys
root = sys.argv[1] if len(sys.argv) > 1 else '.'
cov = open(root+'/data/cycle-three/CubicCovariants.lean').read()
coeffs = {}
for full, base, i, j, k in re.findall(r'def ((H|Theta|J)Coeff(\d)(\d)(\d)) \(a b c d e f : K\)', cov):
    coeffs.setdefault(base, {})[(int(i), int(j), int(k))] = full+'Expr'
# U = a x^3 + b y^3 + c z^3 + d x^2 y + e y^2 z + f z^2 x, coefficients are variables 0..5
coeffs['U'] = {(3,0,0):'(.var 0)', (0,3,0):'(.var 1)', (0,0,3):'(.var 2)',
               (2,1,0):'(.var 3)', (0,2,1):'(.var 4)', (1,0,2):'(.var 5)'}
X, Y, Z = '(.var 6)', '(.var 7)', '(.var 8)'

def term(c, k):
    return c if k == 0 else f'(.mul {c} {Z if k == 1 else f"(.pow {Z} {k})"})'

def horner(base, D):
    """Homogeneous Horner form in y, then x, with explicit z powers; every partial
    sum is weight-homogeneous, as required by the grade checker."""
    C = coeffs[base]; ypart = None
    for j in range(D, -1, -1):
        acc = None
        for i in range(D-j, -1, -1):
            k = D-j-i; c = C.get((i, j, k))
            if acc is None:
                acc = term(c, k) if c else None
            else:
                acc = f'(.mul {acc} {X})'
                if c: acc = f'(.add {acc} {term(c, k)})'
        if ypart is None: ypart = acc
        else:
            ypart = f'(.mul {ypart} {Y})'
            if acc: ypart = f'(.add {ypart} {acc})'
    return ypart

def prod(*fs):
    acc = fs[0]
    for f in fs[1:]: acc = f'(.mul {acc} {f})'
    return acc

def fold(terms):
    acc = terms[0][1]
    for sign, t in terms[1:]:
        acc = f'(.add {acc} {t})' if sign == '+' else f'(.sub {acc} {t})'
    return acc

S, T = 'SExpr', 'TExpr'
outer = [('+', prod('(.num 4)', 'q3')), ('+', prod('(.num 4)', S, 'ThE', 'h4')), ('-', prod(T, 'h6'))]
inner = [('+', prod('(.num 27)', T, 'UE', 'q2')), ('-', prod('(.num 108)', 's3', 'ThE', 'u3')),
         ('+', prod('(.num 54)', S, T, 'ThE', 'u2', 'HE')), ('-', prod('(.num 72)', 's2', 'ThE', 'UE', 'h2')),
         ('-', prod('(.num 18)', T, 'ThE', 'h3')), ('-', prod('(.num 432)', 's4', 'u4', 'HE')),
         ('-', prod('(.num 297)', 's2', T, 'u3', 'h2')), ('-', prod('(.num 108)', 't2', 'u2', 'h3')),
         ('+', prod('(.num 54)', S, T, 'UE', 'h4')), ('-', prod('(.num 16)', 's2', 'h5'))]
hornerBlock = '\n'.join(f'def {name} : Expr :=\n  {horner(base, D)}'
                        for base, name, D in (('U','UE',3), ('H','HE',3), ('Theta','ThE',6), ('J','JE',9)))+'\n\n'
out = 'import Proposal.CubicReflectionData\nimport Proposal.NineKronecker\n\n/-! The cyclic covariant syzygy as one integer identity.\n`U,H,Θ,J` are written homogeneously in `x,y,z` (Horner form), the two sides\nare checked by `CycleNineKronecker.certify`, and the kernel performs the\nbig-number arithmetic. The statement proved is the original `full_syzygy`. -/\nnamespace CycleCovariants.Full\nopen Lean.Grind.CommRing CycleSignedExpressions CycleNineKronecker\nopen CycleCovariants.Reflection\nset_option maxRecDepth 131072\n\n@@HORNER@@def h2 : Expr := .mul HE HE\ndef h3 : Expr := .mul h2 HE\ndef h4 : Expr := .mul h3 HE\ndef h5 : Expr := .mul h4 HE\ndef h6 : Expr := .mul h5 HE\ndef u2 : Expr := .mul UE UE\ndef u3 : Expr := .mul u2 UE\ndef u4 : Expr := .mul u3 UE\ndef q2 : Expr := .mul ThE ThE\ndef q3 : Expr := .mul q2 ThE\ndef s2 : Expr := .mul SExpr SExpr\ndef s3 : Expr := .mul s2 SExpr\ndef s4 : Expr := .mul s3 SExpr\ndef t2 : Expr := .mul TExpr TExpr\ndef LHS : Expr := .mul JE JE\ndef RHS : Expr :=\n  @@RHS@@\ndef base : ℕ := 2^50\n\ntheorem check : certify base LHS RHS=true := by decide +kernel\n\nvariable {K : Type*} [CommRing K]\n\nsection Evaluation\nvariable (a b c d e f x y z : K)\nlemma vars : Lean.RArray.get (context a b c d e f x y z) 0=a ∧ Lean.RArray.get (context a b c d e f x y z) 1=b ∧ Lean.RArray.get (context a b c d e f x y z) 2=c ∧ Lean.RArray.get (context a b c d e f x y z) 3=d ∧ Lean.RArray.get (context a b c d e f x y z) 4=e ∧ Lean.RArray.get (context a b c d e f x y z) 5=f ∧ Lean.RArray.get (context a b c d e f x y z) 6=x ∧ Lean.RArray.get (context a b c d e f x y z) 7=y ∧ Lean.RArray.get (context a b c d e f x y z) 8=z :=\n  ⟨rfl,rfl,rfl,rfl,rfl,rfl,rfl,rfl,rfl⟩\nlemma U_eval : ringEval (fun i => (context a b c d e f x y z).get i) UE=CycleCubic.cubic a b c d e f x y z := by\n  obtain ⟨h0,h1,h2,h3,h4,h5,h6,h7,h8⟩ := vars a b c d e f x y z\n  simp only [UE,ringEval,h0,h1,h2,h3,h4,h5,h6,h7,h8,CycleCubic.cubic]; ring\nlemma H_eval : ringEval (fun i => (context a b c d e f x y z).get i) HE=H a b c d e f x y z := by\n  obtain ⟨-,-,-,-,-,-,h6,h7,h8⟩ := vars a b c d e f x y z\n  simp only [HE,ringEval,h6,h7,h8]\n  simp only [ringEval_denote,@@HDEN@@,H]; ring\nlemma Theta_eval : ringEval (fun i => (context a b c d e f x y z).get i) ThE=Theta a b c d e f x y z := by\n  obtain ⟨-,-,-,-,-,-,h6,h7,h8⟩ := vars a b c d e f x y z\n  simp only [ThE,ringEval,h6,h7,h8]\n  simp only [ringEval_denote,@@THETADEN@@,Theta]; ring\nlemma J_eval : ringEval (fun i => (context a b c d e f x y z).get i) JE=J a b c d e f x y z := by\n  obtain ⟨-,-,-,-,-,-,h6,h7,h8⟩ := vars a b c d e f x y z\n  simp only [JE,ringEval,h6,h7,h8]\n  simp only [ringEval_denote,@@JDEN@@,J]; ring\nend Evaluation\n\ntheorem full_syzygy (a b c d e f x y : K) :\n    let u := CycleCubic.cubic a b c d e f x y 1\n    let h := H a b c d e f x y 1\n    let q := Theta a b c d e f x y 1\n    let j := J a b c d e f x y 1\n    let s := CycleCubic.s4 a b c d e f\n    let t := CycleCubic.s6 a b c d e f\n    j^2=4*q^3+4*s*q*h^4-t*h^6+\n      u*(27*t*u*q^2-108*s^3*q*u^3+54*s*t*q*u^2*h-72*s^2*q*u*h^2-\n        18*t*q*h^3-432*s^4*u^4*h-297*s^2*t*u^3*h^2-108*t^2*u^2*h^3+\n        54*s*t*u*h^4-16*s^2*h^5) := by\n  dsimp only\n  let ctx := context a b c d e f x y 1\n  have hs : ringEval (fun i => ctx.get i) SExpr=CycleCubic.s4 a b c d e f :=\n    (ringEval_denote ctx SExpr).trans (S_denote a b c d e f x y 1)\n  have ht : ringEval (fun i => ctx.get i) TExpr=CycleCubic.s6 a b c d e f :=\n    (ringEval_denote ctx TExpr).trans (T_denote a b c d e f x y 1)\n  have h := certify_sound (fun i => ctx.get i) base LHS RHS check\n  simp only [LHS,RHS,h2,h3,h4,h5,h6,u2,u3,u4,q2,q3,s2,s3,s4,t2,ringEval] at h\n  simp only [ctx,U_eval,H_eval,Theta_eval,J_eval] at h\n  simp only [ctx] at hs ht\n  simp only [hs,ht,Int.cast_ofNat] at h\n  linear_combination h\nend CycleCovariants.Full\n#print axioms CycleCovariants.Full.full_syzygy\n'.replace('@@HORNER@@', hornerBlock)
out = out.replace('@@RHS@@', f'(.add {fold(outer)}\n  (.mul UE {fold(inner)}))\n')
for base, key in (('H','HDEN'), ('Theta','THETADEN'), ('J','JDEN')):
    out = out.replace(f'@@{key}@@', ','.join(v.replace('Expr', '_denote') for v in coeffs[base].values()))
sys.stdout.write(out)
