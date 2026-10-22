import Proposal.CubicReflectionData
import Proposal.SparseSummary

/-! The cyclic covariant syzygy, checked once on integer summaries.
The covariant products are formed by the proved sparse operations, and all
190 coefficients use the encoding certificate of `CubicEncodingCheck`. -/
namespace CycleCovariants.Summary
open Lean.Grind.CommRing CycleSignedExpressions CycleSparseSummary
open CycleCovariants.Reflection
set_option maxRecDepth 131072
set_option maxHeartbeats 4000000

def UE : Poly Expr := [
  ((.var 0),3,0),
  ((.var 3),2,1),
  ((.var 5),1,0),
  ((.var 1),0,3),
  ((.var 4),0,2),
  ((.var 2),0,0)]
def HE : Poly Expr := [
  (HCoeff300Expr,3,0),
  (HCoeff210Expr,2,1),
  (HCoeff201Expr,2,0),
  (HCoeff120Expr,1,2),
  (HCoeff111Expr,1,1),
  (HCoeff102Expr,1,0),
  (HCoeff030Expr,0,3),
  (HCoeff021Expr,0,2),
  (HCoeff012Expr,0,1),
  (HCoeff003Expr,0,0)]
def ThE : Poly Expr := [
  (ThetaCoeff600Expr,6,0),
  (ThetaCoeff510Expr,5,1),
  (ThetaCoeff501Expr,5,0),
  (ThetaCoeff420Expr,4,2),
  (ThetaCoeff411Expr,4,1),
  (ThetaCoeff402Expr,4,0),
  (ThetaCoeff330Expr,3,3),
  (ThetaCoeff321Expr,3,2),
  (ThetaCoeff312Expr,3,1),
  (ThetaCoeff303Expr,3,0),
  (ThetaCoeff240Expr,2,4),
  (ThetaCoeff231Expr,2,3),
  (ThetaCoeff222Expr,2,2),
  (ThetaCoeff213Expr,2,1),
  (ThetaCoeff204Expr,2,0),
  (ThetaCoeff150Expr,1,5),
  (ThetaCoeff141Expr,1,4),
  (ThetaCoeff132Expr,1,3),
  (ThetaCoeff123Expr,1,2),
  (ThetaCoeff114Expr,1,1),
  (ThetaCoeff105Expr,1,0),
  (ThetaCoeff060Expr,0,6),
  (ThetaCoeff051Expr,0,5),
  (ThetaCoeff042Expr,0,4),
  (ThetaCoeff033Expr,0,3),
  (ThetaCoeff024Expr,0,2),
  (ThetaCoeff015Expr,0,1),
  (ThetaCoeff006Expr,0,0)]
def JE : Poly Expr := [
  (JCoeff900Expr,9,0),
  (JCoeff810Expr,8,1),
  (JCoeff801Expr,8,0),
  (JCoeff720Expr,7,2),
  (JCoeff711Expr,7,1),
  (JCoeff702Expr,7,0),
  (JCoeff630Expr,6,3),
  (JCoeff621Expr,6,2),
  (JCoeff612Expr,6,1),
  (JCoeff603Expr,6,0),
  (JCoeff540Expr,5,4),
  (JCoeff531Expr,5,3),
  (JCoeff522Expr,5,2),
  (JCoeff513Expr,5,1),
  (JCoeff504Expr,5,0),
  (JCoeff450Expr,4,5),
  (JCoeff441Expr,4,4),
  (JCoeff432Expr,4,3),
  (JCoeff423Expr,4,2),
  (JCoeff414Expr,4,1),
  (JCoeff405Expr,4,0),
  (JCoeff360Expr,3,6),
  (JCoeff351Expr,3,5),
  (JCoeff342Expr,3,4),
  (JCoeff333Expr,3,3),
  (JCoeff324Expr,3,2),
  (JCoeff315Expr,3,1),
  (JCoeff306Expr,3,0),
  (JCoeff270Expr,2,7),
  (JCoeff261Expr,2,6),
  (JCoeff252Expr,2,5),
  (JCoeff243Expr,2,4),
  (JCoeff234Expr,2,3),
  (JCoeff225Expr,2,2),
  (JCoeff216Expr,2,1),
  (JCoeff207Expr,2,0),
  (JCoeff180Expr,1,8),
  (JCoeff171Expr,1,7),
  (JCoeff162Expr,1,6),
  (JCoeff153Expr,1,5),
  (JCoeff144Expr,1,4),
  (JCoeff135Expr,1,3),
  (JCoeff126Expr,1,2),
  (JCoeff117Expr,1,1),
  (JCoeff108Expr,1,0),
  (JCoeff090Expr,0,9),
  (JCoeff081Expr,0,8),
  (JCoeff072Expr,0,7),
  (JCoeff063Expr,0,6),
  (JCoeff054Expr,0,5),
  (JCoeff045Expr,0,4),
  (JCoeff036Expr,0,3),
  (JCoeff027Expr,0,2),
  (JCoeff018Expr,0,1),
  (JCoeff009Expr,0,0)]

def Hp2 : Poly Expr := pmul .add .mul HE HE
def Hp3 : Poly Expr := pmul .add .mul HE Hp2
def Hp4 : Poly Expr := pmul .add .mul HE Hp3
def Hp5 : Poly Expr := pmul .add .mul HE Hp4
def Hp6 : Poly Expr := pmul .add .mul HE Hp5
def Up2 : Poly Expr := pmul .add .mul UE UE
def Up3 : Poly Expr := pmul .add .mul UE Up2
def Up4 : Poly Expr := pmul .add .mul UE Up3
def Up5 : Poly Expr := pmul .add .mul UE Up4
def Qp2 : Poly Expr := pmul .add .mul ThE ThE
def Qp3 : Poly Expr := pmul .add .mul ThE Qp2
def LeftE : Poly Expr := pmul .add .mul JE JE
def RightE : Poly Expr :=
  pmerge .add 512 (pscale .mul (.num 4) 0 0 Qp3) <|
  pmerge .add 512 (pscale .mul (.mul (.num 4) SExpr) 0 0 (pmul .add .mul ThE Hp4)) <|
  pmerge .add 512 (pscale .mul (.neg TExpr) 0 0 Hp6) <|
  pmerge .add 512 (pscale .mul (.mul (.num 27) TExpr) 0 0 (pmul .add .mul Up2 Qp2)) <|
  pmerge .add 512 (pscale .mul (.mul (.num (-108)) (.pow SExpr 3)) 0 0 (pmul .add .mul ThE Up4)) <|
  pmerge .add 512 (pscale .mul (.mul (.num 54) (.mul SExpr TExpr)) 0 0
    (pmul .add .mul (pmul .add .mul ThE Up3) HE)) <|
  pmerge .add 512 (pscale .mul (.mul (.num (-72)) (.pow SExpr 2)) 0 0
    (pmul .add .mul (pmul .add .mul ThE Up2) Hp2)) <|
  pmerge .add 512 (pscale .mul (.mul (.num (-18)) TExpr) 0 0
    (pmul .add .mul (pmul .add .mul ThE UE) Hp3)) <|
  pmerge .add 512 (pscale .mul (.mul (.num (-432)) (.pow SExpr 4)) 0 0 (pmul .add .mul Up5 HE)) <|
  pmerge .add 512 (pscale .mul (.mul (.num (-297)) (.mul (.pow SExpr 2) TExpr)) 0 0
    (pmul .add .mul Up4 Hp2)) <|
  pmerge .add 512 (pscale .mul (.mul (.num (-108)) (.pow TExpr 2)) 0 0 (pmul .add .mul Up3 Hp3)) <|
  pmerge .add 512 (pscale .mul (.mul (.num 54) (.mul SExpr TExpr)) 0 0 (pmul .add .mul Up2 Hp4))
    (pscale .mul (.mul (.num (-16)) (.pow SExpr 2)) 0 0 (pmul .add .mul UE Hp5))

/-- The single kernel computation, performed on integer summaries. -/
theorem check : pcheck (pmap summ LeftE) (pmap summ RightE)=true := by
  simp only [LeftE,RightE,Hp2,Hp3,Hp4,Hp5,Hp6,Up2,Up3,Up4,Up5,Qp2,Qp3,
    merge_summ,scale_summ,mul_summ]
  decide +kernel

variable {K : Type*} [CommRing K]

theorem identity (v : ℕ → K) (x y : K) :
    pdenote (ringEval v) x y LeftE=pdenote (ringEval v) x y RightE :=
  pcheck_sound v x y LeftE RightE check

theorem dmul (v : ℕ → K) (x y : K) (p q : Poly Expr) :
    pdenote (ringEval v) x y (pmul .add .mul p q)=
      pdenote (ringEval v) x y p*pdenote (ringEval v) x y q :=
  pmul_denote _ _ _ _ _ (ringEval_add v) (ringEval_mul v) p q
theorem dmerge (v : ℕ → K) (x y : K) (n : ℕ) (p q : Poly Expr) :
    pdenote (ringEval v) x y (pmerge .add n p q)=
      pdenote (ringEval v) x y p+pdenote (ringEval v) x y q :=
  pmerge_denote _ _ _ _ (ringEval_add v) n p q
theorem dscale (v : ℕ → K) (x y : K) (a : Expr) (i j : ℕ) (p : Poly Expr) :
    pdenote (ringEval v) x y (pscale .mul a i j p)=
      ringEval v a*x^i*y^j*pdenote (ringEval v) x y p :=
  pscale_denote _ _ _ _ (ringEval_mul v) a i j p

lemma H_evaluation (a b c d e f x y : K) :
    pdenote (ringEval (fun i => (context a b c d e f x y 1).get i)) x y HE=
      H a b c d e f x y 1 := by
  simp only [HE,pdenote,ringEval_denote,HCoeff300_denote,HCoeff210_denote,HCoeff201_denote,HCoeff120_denote,HCoeff111_denote,HCoeff102_denote,HCoeff030_denote,HCoeff021_denote,HCoeff012_denote,HCoeff003_denote,H,pow_zero,pow_one,one_pow,mul_one,add_zero,add_assoc,mul_assoc]
lemma Theta_evaluation (a b c d e f x y : K) :
    pdenote (ringEval (fun i => (context a b c d e f x y 1).get i)) x y ThE=
      Theta a b c d e f x y 1 := by
  simp only [ThE,pdenote,ringEval_denote,ThetaCoeff600_denote,ThetaCoeff510_denote,ThetaCoeff501_denote,ThetaCoeff420_denote,ThetaCoeff411_denote,ThetaCoeff402_denote,ThetaCoeff330_denote,ThetaCoeff321_denote,ThetaCoeff312_denote,ThetaCoeff303_denote,ThetaCoeff240_denote,ThetaCoeff231_denote,ThetaCoeff222_denote,ThetaCoeff213_denote,ThetaCoeff204_denote,ThetaCoeff150_denote,ThetaCoeff141_denote,ThetaCoeff132_denote,ThetaCoeff123_denote,ThetaCoeff114_denote,ThetaCoeff105_denote,ThetaCoeff060_denote,ThetaCoeff051_denote,ThetaCoeff042_denote,ThetaCoeff033_denote,ThetaCoeff024_denote,ThetaCoeff015_denote,ThetaCoeff006_denote,Theta,pow_zero,pow_one,one_pow,mul_one,add_zero,add_assoc,mul_assoc]
lemma J_evaluation (a b c d e f x y : K) :
    pdenote (ringEval (fun i => (context a b c d e f x y 1).get i)) x y JE=
      J a b c d e f x y 1 := by
  simp only [JE,pdenote,ringEval_denote,JCoeff900_denote,JCoeff810_denote,JCoeff801_denote,JCoeff720_denote,JCoeff711_denote,JCoeff702_denote,JCoeff630_denote,JCoeff621_denote,JCoeff612_denote,JCoeff603_denote,JCoeff540_denote,JCoeff531_denote,JCoeff522_denote,JCoeff513_denote,JCoeff504_denote,JCoeff450_denote,JCoeff441_denote,JCoeff432_denote,JCoeff423_denote,JCoeff414_denote,JCoeff405_denote,JCoeff360_denote,JCoeff351_denote,JCoeff342_denote,JCoeff333_denote,JCoeff324_denote,JCoeff315_denote,JCoeff306_denote,JCoeff270_denote,JCoeff261_denote,JCoeff252_denote,JCoeff243_denote,JCoeff234_denote,JCoeff225_denote,JCoeff216_denote,JCoeff207_denote,JCoeff180_denote,JCoeff171_denote,JCoeff162_denote,JCoeff153_denote,JCoeff144_denote,JCoeff135_denote,JCoeff126_denote,JCoeff117_denote,JCoeff108_denote,JCoeff090_denote,JCoeff081_denote,JCoeff072_denote,JCoeff063_denote,JCoeff054_denote,JCoeff045_denote,JCoeff036_denote,JCoeff027_denote,JCoeff018_denote,JCoeff009_denote,J,pow_zero,pow_one,one_pow,mul_one,add_zero,add_assoc,mul_assoc]
lemma U_evaluation (a b c d e f x y : K) :
    pdenote (ringEval (fun i => (context a b c d e f x y 1).get i)) x y UE=
      CycleCubic.cubic a b c d e f x y 1 := by
  have h0 : ringEval (fun i => (context a b c d e f x y 1).get i) (.var 0)=a := rfl
  have h1 : ringEval (fun i => (context a b c d e f x y 1).get i) (.var 1)=b := rfl
  have h2 : ringEval (fun i => (context a b c d e f x y 1).get i) (.var 2)=c := rfl
  have h3 : ringEval (fun i => (context a b c d e f x y 1).get i) (.var 3)=d := rfl
  have h4 : ringEval (fun i => (context a b c d e f x y 1).get i) (.var 4)=e := rfl
  have h5 : ringEval (fun i => (context a b c d e f x y 1).get i) (.var 5)=f := rfl
  simp only [UE,pdenote,
    h0,h1,h2,h3,h4,h5,CycleCubic.cubic]
  ring

theorem full_syzygy (a b c d e f x y : K) :
    let u := CycleCubic.cubic a b c d e f x y 1
    let h := H a b c d e f x y 1
    let q := Theta a b c d e f x y 1
    let j := J a b c d e f x y 1
    let s := CycleCubic.s4 a b c d e f
    let t := CycleCubic.s6 a b c d e f
    j^2=4*q^3+4*s*q*h^4-t*h^6+
      u*(27*t*u*q^2-108*s^3*q*u^3+54*s*t*q*u^2*h-72*s^2*q*u*h^2-
        18*t*q*h^3-432*s^4*u^4*h-297*s^2*t*u^3*h^2-108*t^2*u^2*h^3+
        54*s*t*u*h^4-16*s^2*h^5) := by
  dsimp only
  let ctx := context a b c d e f x y 1
  have hs : ringEval (fun i => ctx.get i) SExpr=CycleCubic.s4 a b c d e f :=
    (ringEval_denote ctx SExpr).trans (S_denote a b c d e f x y 1)
  have ht : ringEval (fun i => ctx.get i) TExpr=CycleCubic.s6 a b c d e f :=
    (ringEval_denote ctx TExpr).trans (T_denote a b c d e f x y 1)
  have h := identity (fun i => ctx.get i) x y
  simp only [LeftE,RightE,Hp2,Hp3,Hp4,Hp5,Hp6,Up2,Up3,Up4,Up5,Qp2,Qp3,
    dmul,dmerge,dscale] at h
  simp only [ctx,H_evaluation,Theta_evaluation,J_evaluation,U_evaluation] at h
  simp only [ctx] at hs ht
  simp only [hs,ht,ringEval,Int.cast_ofNat,Int.cast_neg,pow_zero,mul_one] at h
  linear_combination h
end CycleCovariants.Summary
#print axioms CycleCovariants.Summary.full_syzygy
