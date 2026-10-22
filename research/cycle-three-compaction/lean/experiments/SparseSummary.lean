import Proposal.CubicEncodingCheck

/-! Sparse bivariate polynomials over an arbitrary coefficient type, and an
integer summary of ring expressions that commutes with their operations.

Products are formed once on summaries (weight grade, positive and negative
mass at one, value at the Kronecker point), so the kernel never evaluates
expanded product expressions. Soundness reduces every coefficient to the
encoding certificate `identity_of_encoding_context`. -/
namespace CycleSparseSummary
open Lean.Grind.CommRing CycleSignedExpressions
open CycleCubicEncodingCheck CycleCubicGradeCheck CycleCubicKronecker

abbrev Poly (α : Type) := List (α × ℕ × ℕ)

section Generic
variable {α β : Type}

def pmerge (add : α → α → α) : ℕ → Poly α → Poly α → Poly α
  | 0,p,q => p++q
  | _+1,[],q => q
  | _+1,p,[] => p
  | n+1,(a,i,j)::p,(b,k,l)::q =>
    if i=k ∧ j=l then (add a b,i,j)::pmerge add n p q
    else if 32*i+j>32*k+l then (a,i,j)::pmerge add n p ((b,k,l)::q)
    else (b,k,l)::pmerge add n ((a,i,j)::p) q

def pscale (mul : α → α → α) (a : α) (i j : ℕ) : Poly α → Poly α
  | [] => []
  | (b,k,l)::q => (mul a b,i+k,j+l)::pscale mul a i j q

def pmul (add mul : α → α → α) : Poly α → Poly α → Poly α
  | [],_ => []
  | (a,i,j)::p,q => pmerge add 512 (pscale mul a i j q) (pmul add mul p q)

def pmap (f : α → β) (p : Poly α) : Poly β := p.map fun t => (f t.1,t.2.1,t.2.2)

section Naturality
variable (add mul : α → α → α) (add' mul' : β → β → β) (f : α → β)

theorem pmerge_map (hadd : ∀ a b,f (add a b)=add' (f a) (f b)) :
    ∀ n p q,pmap f (pmerge add n p q)=pmerge add' n (pmap f p) (pmap f q)
  | 0,p,q => by simp [pmerge,pmap]
  | _+1,[],q => by simp [pmerge,pmap]
  | _+1,_::_,[] => by simp [pmerge,pmap]
  | n+1,(a,i,j)::p,(b,k,l)::q => by
    have ih := pmerge_map hadd n
    have hc : ∀ (c : α) (r s : ℕ) (t : Poly α),pmap f ((c,r,s)::t)=(f c,r,s)::pmap f t :=
      fun _ _ _ _ => rfl
    rw [hc,hc]
    simp only [pmerge]
    split_ifs
    · rw [hc,ih,hadd]
    · rw [hc,ih,hc]
    · rw [hc,ih,hc]

theorem pscale_map (hmul : ∀ a b,f (mul a b)=mul' (f a) (f b)) (a : α) (i j : ℕ) :
    ∀ p,pmap f (pscale mul a i j p)=pscale mul' (f a) i j (pmap f p)
  | [] => rfl
  | (b,k,l)::q => by
    have := pscale_map hmul a i j q
    simp only [pscale,pmap,List.map_cons] at this ⊢
    simp [this,hmul]

theorem pmul_map (hadd : ∀ a b,f (add a b)=add' (f a) (f b))
    (hmul : ∀ a b,f (mul a b)=mul' (f a) (f b)) :
    ∀ p q,pmap f (pmul add mul p q)=pmul add' mul' (pmap f p) (pmap f q)
  | [],q => rfl
  | (a,i,j)::p,q => by
    have := pmul_map hadd hmul p q
    simp only [pmul]
    rw [pmerge_map add add' f hadd,pscale_map mul mul' f hmul,this]
    rfl
end Naturality

section Denote
variable {K : Type*} [CommRing K] (add mul : α → α → α) (ev : α → K) (x y : K)

def pdenote : Poly α → K
  | [] => 0
  | (a,i,j)::p => ev a*x^i*y^j+pdenote p

theorem pmerge_denote (hadd : ∀ a b,ev (add a b)=ev a+ev b) :
    ∀ n p q,pdenote ev x y (pmerge add n p q)=pdenote ev x y p+pdenote ev x y q
  | 0,p,q => by
    induction p with
    | nil => simp [pmerge,pdenote]
    | cons t p ih =>
      obtain ⟨a,i,j⟩ := t
      simp only [pmerge,List.cons_append,pdenote] at ih ⊢
      rw [ih]; ring
  | _+1,[],q => by simp [pmerge,pdenote]
  | _+1,_::_,[] => by simp [pmerge,pdenote]
  | n+1,(a,i,j)::p,(b,k,l)::q => by
    have ih := pmerge_denote hadd n
    simp only [pmerge]
    split_ifs with h h'
    · obtain ⟨rfl,rfl⟩ := h
      simp only [pdenote,ih,hadd]; ring
    · simp only [pdenote,ih]; ring
    · simp only [pdenote,ih]; ring

theorem pscale_denote (hmul : ∀ a b,ev (mul a b)=ev a*ev b) (a : α) (i j : ℕ) :
    ∀ p,pdenote ev x y (pscale mul a i j p)=ev a*x^i*y^j*pdenote ev x y p
  | [] => by simp [pscale,pdenote]
  | (b,k,l)::q => by
    simp only [pscale,pdenote,pscale_denote hmul a i j q,hmul,pow_add]; ring

theorem pmul_denote (hadd : ∀ a b,ev (add a b)=ev a+ev b)
    (hmul : ∀ a b,ev (mul a b)=ev a*ev b) :
    ∀ p q,pdenote ev x y (pmul add mul p q)=pdenote ev x y p*pdenote ev x y q
  | [],q => by simp [pmul,pdenote]
  | (a,i,j)::p,q => by
    simp only [pmul]
    rw [pmerge_denote add ev x y hadd,pscale_denote mul ev x y hmul,
      pmul_denote hadd hmul p q]
    simp only [pdenote]; ring
end Denote
end Generic

/-! Integer summaries of ring expressions. -/
structure Summ where
  grade : Option Grade
  pos : ℕ
  neg : ℕ
  val : ℤ

def base : ℕ := 2^42
def sadd (s t : Summ) : Summ := ⟨gradeAdd s.grade t.grade,s.pos+t.pos,s.neg+t.neg,s.val+t.val⟩
def smul (s t : Summ) : Summ :=
  ⟨gradeMul s.grade t.grade,s.pos*t.pos+s.neg*t.neg,s.pos*t.neg+s.neg*t.pos,s.val*t.val⟩
def summ (e : Expr) : Summ :=
  ⟨grade e,(splitEval (value 1) e).1,(splitEval (value 1) e).2,
    ringEval (fun i => (value base i : ℤ)) e⟩

theorem summ_add (a b : Expr) : summ (.add a b)=sadd (summ a) (summ b) := rfl
theorem summ_mul (a b : Expr) : summ (.mul a b)=smul (summ a) (summ b) := rfl

def coeffOK (s t : Summ) : Bool :=
  match s.grade with
  | some w =>
    decide (t.grade=some w) && decide (w.1<38 ∧ w.2.1<38 ∧ w.2.2<38) &&
      decide (s.pos+t.neg<base) && decide (t.pos+s.neg<base) && decide (s.val=t.val)
  | none => false

def pcheck : Poly Summ → Poly Summ → Bool
  | [],[] => true
  | (s,i,j)::p,(t,k,l)::q => decide (i=k) && decide (j=l) && coeffOK s t && pcheck p q
  | _,_ => false

variable {K : Type*} [CommRing K]

theorem coeffOK_sound (v : ℕ → K) (e f : Expr) (h : coeffOK (summ e) (summ f)=true) :
    ringEval v e=ringEval v f := by
  unfold coeffOK at h
  split at h
  · rename_i w hw
    simp only [Bool.and_eq_true,decide_eq_true_eq] at h
    obtain ⟨⟨⟨⟨hf,hb⟩,hp⟩,hn⟩,hnum⟩ := h
    exact identity_of_encoding_context v e f w base hb (by decide) hw hf hp hn hnum
  · exact absurd h (by simp)

theorem pcheck_sound (v : ℕ → K) (x y : K) :
    ∀ p q : Poly Expr,pcheck (pmap summ p) (pmap summ q)=true →
      pdenote (ringEval v) x y p=pdenote (ringEval v) x y q
  | [],[],_ => rfl
  | (a,i,j)::p,(b,k,l)::q,h => by
    simp only [pmap,List.map_cons,pcheck,Bool.and_eq_true,decide_eq_true_eq] at h
    obtain ⟨⟨⟨rfl,rfl⟩,hab⟩,hpq⟩ := h
    simp only [pdenote,coeffOK_sound v a b hab,pcheck_sound v x y p q hpq]
  | [],_::_,h => by simp [pmap,pcheck] at h
  | _::_,[],h => by simp [pmap,pcheck] at h

/-! Specializations used to move a whole product program onto summaries. -/
theorem merge_summ (n : ℕ) (p q : Poly Expr) :
    pmap summ (pmerge .add n p q)=pmerge sadd n (pmap summ p) (pmap summ q) :=
  pmerge_map _ _ _ summ_add n p q
theorem scale_summ (a : Expr) (i j : ℕ) (p : Poly Expr) :
    pmap summ (pscale .mul a i j p)=pscale smul (summ a) i j (pmap summ p) :=
  pscale_map _ _ _ summ_mul a i j p
theorem mul_summ (p q : Poly Expr) :
    pmap summ (pmul .add .mul p q)=pmul sadd smul (pmap summ p) (pmap summ q) :=
  pmul_map _ _ _ _ _ summ_add summ_mul p q

theorem ringEval_add (v : ℕ → K) (a b : Expr) :
    ringEval v (.add a b)=ringEval v a+ringEval v b := rfl
theorem ringEval_mul (v : ℕ → K) (a b : Expr) :
    ringEval v (.mul a b)=ringEval v a*ringEval v b := rfl
end CycleSparseSummary
#print axioms CycleSparseSummary.pcheck_sound
