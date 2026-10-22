import Proposal.SignedExpressions
import Proposal.PolynomialDigitCheck
import Mathlib.RingTheory.MvPolynomial.WeightedHomogeneous

/-! A single positional encoding for whole covariant identities.

Variables are the cubic coefficients `a,b,c,d,e,f` (0–5) and the homogeneous
coordinates `x,y,z` (6–8). The four weights below make the equation and every
covariant of the cyclic cubic homogeneous. At fixed weight with all four
components bounded, the exponents of `a,b,c,z` are determined by those of
`d,e,f,x,y`, which are at most 18. Substituting `a=b=c=z=1`,
`d,e,f,x,y = X, X^19, X^361, X^6859, X^130321` is therefore injective on that
component, and evaluation at a base exceeding the coefficient mass is
injective by `CyclePolynomialDigits`. A whole identity is then one integer
equality, computed by the kernel's big-number arithmetic. -/
noncomputable section
namespace CycleNineKronecker
open MvPolynomial Lean.Grind.CommRing CycleSignedExpressions
open scoped BigOperators

abbrev Grade := ℕ × ℕ × ℕ × ℕ
abbrev Exponent := Fin 9 →₀ ℕ
def weights : Fin 9 → Grade :=
  ![(3,0,0,0),(0,3,0,0),(0,0,3,0),(2,1,0,0),(0,2,1,0),(1,0,2,0),
    (0,1,1,1),(1,0,1,1),(1,1,0,1)]
def strides : Fin 9 → ℕ := ![0,0,0,1,19,361,6859,130321,0]
def code (m : Exponent) : ℕ := ∑ i,strides i*m i

lemma code_eq (m : Exponent) : code m=m 3+19*m 4+361*m 5+6859*m 6+130321*m 7 := by
  simp [code,Fin.sum_univ_succ,strides,add_assoc]

lemma weight_formula (m : Exponent) :
    Finsupp.weight weights m=
      (3*m 0+2*m 3+m 5+m 7+m 8,3*m 1+m 3+2*m 4+m 6+m 8,
        3*m 2+m 4+2*m 5+m 6+m 7,m 6+m 7+m 8) := by
  rw [Finsupp.weight_apply,Finsupp.sum_fintype]
  · simp only [weights,Fin.sum_univ_succ,Fin.sum_univ_zero,Matrix.cons_val_zero,
      Matrix.cons_val_succ,Prod.smul_mk,Prod.mk_add_mk,smul_eq_mul,add_zero,mul_zero,
      mul_one,Fin.succ_zero_eq_one,Fin.succ_one_eq_two]
    simp only [Prod.mk.injEq]
    refine ⟨?_,?_,?_,?_⟩ <;> simp <;> ring
  · intro; simp

def Bounded (w : Grade) : Prop := w.1<38 ∧ w.2.1<38 ∧ w.2.2.1<38 ∧ w.2.2.2<19
instance (w : Grade) : Decidable (Bounded w) := by unfold Bounded; infer_instance

lemma digit_bound (a b c d e f x y z u v w t : ℕ)
    (h1 : 3*a+2*d+f+y+z=u) (h2 : 3*b+d+2*e+x+z=v) (h3 : 3*c+e+2*f+x+y=w)
    (h4 : x+y+z=t) (hu : u<38) (hv : v<38) (hw : w<38) (ht : t<19) :
    d<19 ∧ e<19 ∧ f<19 ∧ x<19 := by omega

lemma digits_unique (d e f x y d' e' f' x' y' : ℕ)
    (hd : d<19) (he : e<19) (hf : f<19) (hx : x<19)
    (hd' : d'<19) (he' : e'<19) (hf' : f'<19) (hx' : x'<19)
    (h : d+19*e+361*f+6859*x+130321*y=d'+19*e'+361*f'+6859*x'+130321*y') :
    d=d' ∧ e=e' ∧ f=f' ∧ x=x' ∧ y=y' := by omega

lemma recover (a b c d e f x y z a' b' c' d' e' f' x' y' z' u v w t : ℕ)
    (h1 : 3*a+2*d+f+y+z=u) (h2 : 3*b+d+2*e+x+z=v) (h3 : 3*c+e+2*f+x+y=w) (h4 : x+y+z=t)
    (h1' : 3*a'+2*d'+f'+y'+z'=u) (h2' : 3*b'+d'+2*e'+x'+z'=v) (h3' : 3*c'+e'+2*f'+x'+y'=w)
    (h4' : x'+y'+z'=t) (hd : d=d') (he : e=e') (hf : f=f') (hx : x=x') (hy : y=y') :
    a=a' ∧ b=b' ∧ c=c' ∧ z=z' := by subst hd he hf hx hy; omega

lemma code_injective_on_grade (w : Grade) (hw : Bounded w) :
    Set.InjOn code {m : Exponent | Finsupp.weight weights m=w} := by
  intro m hm n hn he
  rw [Set.mem_setOf_eq,weight_formula] at hm hn
  obtain ⟨w1,w2,w3,w4⟩ := w
  simp only [Prod.mk.injEq] at hm hn
  obtain ⟨hm1,hm2,hm3,hm4⟩ := hm
  obtain ⟨hn1,hn2,hn3,hn4⟩ := hn
  obtain ⟨b1,b2,b3,b4⟩ := hw
  simp only at b1 b2 b3 b4
  rw [code_eq,code_eq] at he
  obtain ⟨d1,d2,d3,d4⟩ := digit_bound _ _ _ _ _ _ _ _ _ _ _ _ _ hm1 hm2 hm3 hm4 b1 b2 b3 b4
  obtain ⟨e1,e2,e3,e4⟩ := digit_bound _ _ _ _ _ _ _ _ _ _ _ _ _ hn1 hn2 hn3 hn4 b1 b2 b3 b4
  obtain ⟨h3,h4,h5,h6,h7⟩ := digits_unique _ _ _ _ _ _ _ _ _ _ d1 d2 d3 d4 e1 e2 e3 e4 he
  obtain ⟨h0,h1,h2,h8⟩ := recover _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    hm1 hm2 hm3 hm4 hn1 hn2 hn3 hn4 h3 h4 h5 h6 h7
  ext i
  fin_cases i <;> assumption

variable {R : Type*} [CommSemiring R]
def encode : MvPolynomial (Fin 9) R →+* Polynomial R :=
  MvPolynomial.eval₂Hom Polynomial.C fun i => Polynomial.X^strides i

lemma encode_monomial (m : Exponent) (r : R) :
    encode (MvPolynomial.monomial m r)=Polynomial.monomial (code m) r := by
  rw [encode,MvPolynomial.eval₂Hom_monomial,
    Finsupp.prod_fintype _ _ (fun i => pow_zero _),
    Finset.prod_congr rfl (fun i _ => (pow_mul Polynomial.X (strides i) (m i)).symm),
    Finset.prod_pow_eq_pow_sum,Polynomial.C_mul_X_pow_eq_monomial]
  rfl

lemma encode_finsupp (p : MvPolynomial (Fin 9) R) :
    (encode p).toFinsupp.coeff=Finsupp.mapDomain code (AddMonoidAlgebra.coeff p) := by
  induction p using MvPolynomial.induction_on' with
  | monomial m r =>
    rw [encode_monomial]
    simp [Polynomial.toFinsupp_monomial,MvPolynomial.monomial,Finsupp.mapDomain_single]
  | add p q hp hq =>
    simp only [map_add,Polynomial.toFinsupp_add,AddMonoidAlgebra.coeff_add,
      Finsupp.mapDomain_add,hp,hq]

theorem encode_injective_on_grade (w : Grade) (hw : Bounded w)
    (p q : MvPolynomial (Fin 9) R)
    (hp : p.IsWeightedHomogeneous weights w) (hq : q.IsWeightedHomogeneous weights w)
    (he : encode p=encode q) : p=q := by
  classical
  let S : Set Exponent := {m | Finsupp.weight weights m=w}
  have hsP : (p.support : Set Exponent)⊆S := fun m hm => hp (mem_support_iff.mp hm)
  have hsQ : (q.support : Set Exponent)⊆S := fun m hm => hq (mem_support_iff.mp hm)
  have heF := congrArg (fun r : Polynomial R => r.toFinsupp.coeff) he
  rw [encode_finsupp,encode_finsupp] at heF
  ext m
  by_cases hm : m∈S
  · have h := congrArg (fun f : ℕ →₀ R => f (code m)) heF
    rw [Finsupp.mapDomain_apply' S (AddMonoidAlgebra.coeff p) hsP (code_injective_on_grade w hw) hm,
      Finsupp.mapDomain_apply' S (AddMonoidAlgebra.coeff q) hsQ (code_injective_on_grade w hw) hm] at h
    exact h
  · have hp0 : coeff m p=0 := by
      by_contra h; exact hm (hp h)
    have hq0 : coeff m q=0 := by
      by_contra h; exact hm (hq h)
    rw [hp0,hq0]

/-! A proved syntactic weight checker. -/
def gradeAdd : Option Grade → Option Grade → Option Grade
  | some a,some b => if a=b then some a else none
  | _,_ => none
def gradeMul : Option Grade → Option Grade → Option Grade
  | some a,some b => some (a+b)
  | _,_ => none
def grade : Expr → Option Grade
  | .num _ | .natCast _ | .intCast _ => some 0
  | .var i => if h : i<9 then some (weights ⟨i,h⟩) else none
  | .neg a => grade a
  | .add a b | .sub a b => gradeAdd (grade a) (grade b)
  | .mul a b => gradeMul (grade a) (grade b)
  | .pow a n => (grade a).map (n • ·)

def polyVariable (i : ℕ) : MvPolynomial (Fin 9) ℕ :=
  if h : i<9 then X ⟨i,h⟩ else 0

def Homogeneous (p : MvPolynomial (Fin 9) ℕ × MvPolynomial (Fin 9) ℕ)
    (w : Grade) : Prop :=
  p.1.IsWeightedHomogeneous weights w ∧ p.2.IsWeightedHomogeneous weights w

lemma homogeneous_number (n : Int) : Homogeneous (numberPair n) 0 := by
  cases n with
  | ofNat n =>
    exact ⟨by simpa [numberPair] using isWeightedHomogeneous_C weights n,
      isWeightedHomogeneous_zero ℕ weights 0⟩
  | negSucc n =>
    exact ⟨isWeightedHomogeneous_zero ℕ weights 0,
      by simpa [numberPair] using isWeightedHomogeneous_C weights (n+1)⟩

lemma homogeneous_mul {p q : MvPolynomial (Fin 9) ℕ × MvPolynomial (Fin 9) ℕ}
    {u v : Grade} (hp : Homogeneous p u) (hq : Homogeneous q v) :
    Homogeneous (pairMul p q) (u+v) :=
  ⟨(hp.1.mul hq.1).add (hp.2.mul hq.2),
    (hp.1.mul hq.2).add (hp.2.mul hq.1)⟩

lemma homogeneous_pow {p : MvPolynomial (Fin 9) ℕ × MvPolynomial (Fin 9) ℕ}
    {u : Grade} (hp : Homogeneous p u) (n : ℕ) :
    Homogeneous (pairPow p n) (n • u) := by
  induction n with
  | zero =>
    simpa only [Homogeneous,pairPow,zero_nsmul] using
      And.intro (isWeightedHomogeneous_one ℕ weights) (isWeightedHomogeneous_zero ℕ weights 0)
  | succ n ih =>
    simpa only [pairPow,succ_nsmul',add_comm] using homogeneous_mul hp ih

theorem grade_sound (e : Expr) (w : Grade) (h : grade e=some w) :
    Homogeneous (splitEval polyVariable e) w := by
  induction e generalizing w with
  | num n => cases Option.some.inj h; exact homogeneous_number n
  | intCast n => cases Option.some.inj h; exact homogeneous_number n
  | natCast n =>
    cases Option.some.inj h
    exact ⟨by simpa [splitEval] using isWeightedHomogeneous_C weights n,
      isWeightedHomogeneous_zero ℕ weights 0⟩
  | var i =>
    simp only [grade] at h
    split_ifs at h with hi
    · cases Option.some.inj h
      exact ⟨by simpa [splitEval,polyVariable,hi] using
        isWeightedHomogeneous_X (R:=ℕ) weights (⟨i,hi⟩ : Fin 9),
        isWeightedHomogeneous_zero ℕ weights _⟩
  | neg a ih => exact ⟨(ih w h).2,(ih w h).1⟩
  | add a b ha hb =>
    cases hga : grade a with
    | none => simp [grade,gradeAdd,hga] at h
    | some ga =>
      cases hgb : grade b with
      | none => simp [grade,gradeAdd,hga,hgb] at h
      | some gb =>
        simp only [grade,hga,hgb,gradeAdd] at h
        split_ifs at h with he
        · subst gb
          have hw : ga=w := Option.some.inj h
          subst w
          exact ⟨((ha ga hga).1).add ((hb ga hgb).1),
            ((ha ga hga).2).add ((hb ga hgb).2)⟩
  | sub a b ha hb =>
    cases hga : grade a with
    | none => simp [grade,gradeAdd,hga] at h
    | some ga =>
      cases hgb : grade b with
      | none => simp [grade,gradeAdd,hga,hgb] at h
      | some gb =>
        simp only [grade,hga,hgb,gradeAdd] at h
        split_ifs at h with he
        · subst gb
          have hw : ga=w := Option.some.inj h
          subst w
          exact ⟨((ha ga hga).1).add ((hb ga hgb).2),
            ((ha ga hga).2).add ((hb ga hgb).1)⟩
  | mul a b ha hb =>
    cases hga : grade a with
    | none => simp [grade,gradeMul,hga] at h
    | some ga =>
      cases hgb : grade b with
      | none => simp [grade,gradeMul,hga,hgb] at h
      | some gb =>
        simp only [grade,hga,hgb,gradeMul] at h
        have hw : ga+gb=w := Option.some.inj h
        subst w
        exact homogeneous_mul (ha ga hga) (hb gb hgb)
  | pow a n ha =>
    cases hga : grade a with
    | none => simp [grade,hga] at h
    | some ga =>
      simp only [grade,hga,Option.map_some] at h
      have hw : n • ga=w := Option.some.inj h
      subst w
      exact homogeneous_pow (ha ga hga) n

lemma gradeAdd_valid {a b : Option Grade} (h : gradeAdd a b≠none) :
    a≠none ∧ b≠none := by
  cases a <;> cases b <;> simp_all [gradeAdd]
lemma gradeMul_valid {a b : Option Grade} (h : gradeMul a b≠none) :
    a≠none ∧ b≠none := by
  cases a <;> cases b <;> simp_all [gradeMul]

theorem ringEval_congr_of_grade {K : Type*} [CommRing K]
    (e : Expr) (v v' : ℕ → K) (hv : ∀ i,i<9 → v i=v' i)
    (h : grade e≠none) : ringEval v e=ringEval v' e := by
  induction e with
  | num n => rfl
  | natCast n => rfl
  | intCast n => rfl
  | var i =>
    by_cases hi : i<9
    · exact hv i hi
    · simp [grade,hi] at h
  | neg a ha => exact congrArg Neg.neg (ha h)
  | add a b ha hb =>
    obtain ⟨h₁,h₂⟩ := gradeAdd_valid h
    exact congrArg₂ HAdd.hAdd (ha h₁) (hb h₂)
  | sub a b ha hb =>
    obtain ⟨h₁,h₂⟩ := gradeAdd_valid h
    exact congrArg₂ HSub.hSub (ha h₁) (hb h₂)
  | mul a b ha hb =>
    obtain ⟨h₁,h₂⟩ := gradeMul_valid h
    exact congrArg₂ HMul.hMul (ha h₁) (hb h₂)
  | pow a n ha =>
    have h' : grade a≠none := by
      intro hz
      simp [grade,hz] at h
    exact congrArg (fun z : K => z^n) (ha h')

/-! The encoding certificate. -/
def value (b : ℕ) (i : ℕ) : ℕ := if h : i<9 then b^strides ⟨i,h⟩ else 0
def evaluation (b : ℕ) : MvPolynomial (Fin 9) ℕ →+* ℕ :=
  (Polynomial.evalRingHom b).comp encode

lemma evaluation_variable (b i : ℕ) : evaluation b (polyVariable i)=value b i := by
  unfold polyVariable value
  split_ifs with hi
  · simp [evaluation,encode]
  · simp

lemma evaluation_parts (b : ℕ) (e : Expr) :
    (evaluation b (splitEval polyVariable e).1,evaluation b (splitEval polyVariable e).2)=
      splitEval (value b) e := by
  simpa only [evaluation_variable] using splitEval_map (evaluation b) polyVariable e

lemma cast_parts (b : ℕ) (e : Expr) :
    ((splitEval (value b) e).1 : ℤ)-(splitEval (value b) e).2=
      ringEval (fun i => (value b i : ℤ)) e := by
  have h' := splitEval_map (Nat.castRingHom ℤ) (value b) e
  change (((splitEval (value b) e).1 : ℤ),((splitEval (value b) e).2 : ℤ))=
    splitEval (fun i => (value b i : ℤ)) e at h'
  have hd := splitEval_difference (fun i => (value b i : ℤ)) e
  rw [← h'] at hd
  exact hd

def positiveMass (e f : Expr) : ℕ :=
  (splitEval (value 1) e).1+(splitEval (value 1) f).2
def negativeMass (e f : Expr) : ℕ :=
  (splitEval (value 1) f).1+(splitEval (value 1) e).2

theorem encoded_identity (e f : Expr) (w : Grade) (b : ℕ)
    (hw : Bounded w) (hb : 0<b)
    (he : grade e=some w) (hf : grade f=some w)
    (hp : positiveMass e f<b) (hn : negativeMass e f<b)
    (hnum : ringEval (fun i => (value b i : ℤ)) e=
      ringEval (fun i => (value b i : ℤ)) f) :
    (splitEval polyVariable e).1+(splitEval polyVariable f).2=
      (splitEval polyVariable f).1+(splitEval polyVariable e).2 := by
  let p := (splitEval polyVariable e).1+(splitEval polyVariable f).2
  let q := (splitEval polyVariable f).1+(splitEval polyVariable e).2
  have hew := grade_sound e w he
  have hfw := grade_sound f w hf
  apply encode_injective_on_grade w hw p q (hew.1.add hfw.2) (hfw.1.add hew.2)
  apply CyclePolynomialDigits.eval_injective_bounded b hb
  · change evaluation 1 p<b
    have he1 := congrArg Prod.fst (evaluation_parts 1 e)
    have hf2 := congrArg Prod.snd (evaluation_parts 1 f)
    dsimp only at he1 hf2
    simpa only [p,map_add,he1,hf2,positiveMass] using hp
  · change evaluation 1 q<b
    have hf1 := congrArg Prod.fst (evaluation_parts 1 f)
    have he2 := congrArg Prod.snd (evaluation_parts 1 e)
    dsimp only at hf1 he2
    simpa only [q,map_add,hf1,he2,negativeMass] using hn
  · change evaluation b p=evaluation b q
    have he1 := congrArg Prod.fst (evaluation_parts b e)
    have he2 := congrArg Prod.snd (evaluation_parts b e)
    have hf1 := congrArg Prod.fst (evaluation_parts b f)
    have hf2 := congrArg Prod.snd (evaluation_parts b f)
    dsimp only at he1 he2 hf1 hf2
    simp only [p,q,map_add,he1,he2,hf1,hf2]
    rw [← cast_parts,← cast_parts] at hnum
    omega

theorem identity_of_encoding {K : Type*} [CommRing K]
    (v : Fin 9 → K) (e f : Expr) (w : Grade) (b : ℕ)
    (hw : Bounded w) (hb : 0<b)
    (he : grade e=some w) (hf : grade f=some w)
    (hp : positiveMass e f<b) (hn : negativeMass e f<b)
    (hnum : ringEval (fun i => (value b i : ℤ)) e=
      ringEval (fun i => (value b i : ℤ)) f) :
    ringEval (fun i => if h : i<9 then v ⟨i,h⟩ else 0) e=
      ringEval (fun i => if h : i<9 then v ⟨i,h⟩ else 0) f := by
  let φ : MvPolynomial (Fin 9) ℕ →+* K := eval₂Hom (Nat.castRingHom K) v
  have hvar (i : ℕ) : φ (polyVariable i)=(if h : i<9 then v ⟨i,h⟩ else 0) := by
    by_cases hi : i<9 <;> simp [φ,polyVariable,hi]
  have h := congrArg φ (encoded_identity e f w b hw hb he hf hp hn hnum)
  rw [map_add,map_add] at h
  have he' := splitEval_map φ polyVariable e
  have hf' := splitEval_map φ polyVariable f
  have heDiff := splitEval_difference (fun i => φ (polyVariable i)) e
  have hfDiff := splitEval_difference (fun i => φ (polyVariable i)) f
  rw [← he'] at heDiff
  rw [← hf'] at hfDiff
  simp only [hvar] at heDiff hfDiff
  rw [← heDiff,← hfDiff]
  linear_combination h

/-- The complete certificate as a Boolean, for one kernel evaluation. -/
def certify (b : ℕ) (e f : Expr) : Bool :=
  match grade e with
  | some w => decide (grade f=some w) && decide (Bounded w) && decide (0<b) &&
      decide (positiveMass e f<b) && decide (negativeMass e f<b) &&
      decide (ringEval (fun i => (value b i : ℤ)) e=ringEval (fun i => (value b i : ℤ)) f)
  | none => false

theorem certify_sound {K : Type*} [CommRing K] (v : ℕ → K) (b : ℕ) (e f : Expr)
    (h : certify b e f=true) : ringEval v e=ringEval v f := by
  unfold certify at h
  split at h
  · rename_i w he
    simp only [Bool.and_eq_true,decide_eq_true_eq] at h
    obtain ⟨⟨⟨⟨⟨hf,hw⟩,hb⟩,hp⟩,hn⟩,hnum⟩ := h
    let v' : ℕ → K := fun i => if i<9 then v i else 0
    have hv (i : ℕ) (hi : i<9) : v i=v' i := by simp [v',hi]
    have he' : grade e≠none := by rw [he]; exact Option.some_ne_none w
    have hf' : grade f≠none := by rw [hf]; exact Option.some_ne_none w
    rw [ringEval_congr_of_grade e v v' hv he',ringEval_congr_of_grade f v v' hv hf']
    have := identity_of_encoding (fun i : Fin 9 => v i) e f w b hw hb he hf hp hn hnum
    simpa [v'] using this
  · exact absurd h (by simp)
end CycleNineKronecker
#print axioms CycleNineKronecker.certify_sound
