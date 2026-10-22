import Proposal.CubicReflectionData
import Proposal.NineKronecker

/-! The cyclic covariant syzygy as one integer identity.
`U,H,Θ,J` are written homogeneously in `x,y,z` (Horner form), the two sides
are checked by `CycleNineKronecker.certify`, and the kernel performs the
big-number arithmetic. The statement proved is the original `full_syzygy`. -/
namespace CycleCovariants.Full
open Lean.Grind.CommRing CycleSignedExpressions CycleNineKronecker
open CycleCovariants.Reflection
set_option maxRecDepth 131072

def UE : Expr :=
  (.add (.mul (.add (.mul (.add (.mul (.var 1) (.var 7)) (.mul (.var 4) (.var 8))) (.var 7)) (.mul (.mul (.var 3) (.var 6)) (.var 6))) (.var 7)) (.add (.mul (.add (.mul (.mul (.var 0) (.var 6)) (.var 6)) (.mul (.var 5) (.pow (.var 8) 2))) (.var 6)) (.mul (.var 2) (.pow (.var 8) 3))))
def HE : Expr :=
  (.add (.mul (.add (.mul (.add (.mul HCoeff030Expr (.var 7)) (.add (.mul HCoeff120Expr (.var 6)) (.mul HCoeff021Expr (.var 8)))) (.var 7)) (.add (.mul (.add (.mul HCoeff210Expr (.var 6)) (.mul HCoeff111Expr (.var 8))) (.var 6)) (.mul HCoeff012Expr (.pow (.var 8) 2)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul HCoeff300Expr (.var 6)) (.mul HCoeff201Expr (.var 8))) (.var 6)) (.mul HCoeff102Expr (.pow (.var 8) 2))) (.var 6)) (.mul HCoeff003Expr (.pow (.var 8) 3))))
def ThE : Expr :=
  (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul ThetaCoeff060Expr (.var 7)) (.add (.mul ThetaCoeff150Expr (.var 6)) (.mul ThetaCoeff051Expr (.var 8)))) (.var 7)) (.add (.mul (.add (.mul ThetaCoeff240Expr (.var 6)) (.mul ThetaCoeff141Expr (.var 8))) (.var 6)) (.mul ThetaCoeff042Expr (.pow (.var 8) 2)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul ThetaCoeff330Expr (.var 6)) (.mul ThetaCoeff231Expr (.var 8))) (.var 6)) (.mul ThetaCoeff132Expr (.pow (.var 8) 2))) (.var 6)) (.mul ThetaCoeff033Expr (.pow (.var 8) 3)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul (.add (.mul ThetaCoeff420Expr (.var 6)) (.mul ThetaCoeff321Expr (.var 8))) (.var 6)) (.mul ThetaCoeff222Expr (.pow (.var 8) 2))) (.var 6)) (.mul ThetaCoeff123Expr (.pow (.var 8) 3))) (.var 6)) (.mul ThetaCoeff024Expr (.pow (.var 8) 4)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul ThetaCoeff510Expr (.var 6)) (.mul ThetaCoeff411Expr (.var 8))) (.var 6)) (.mul ThetaCoeff312Expr (.pow (.var 8) 2))) (.var 6)) (.mul ThetaCoeff213Expr (.pow (.var 8) 3))) (.var 6)) (.mul ThetaCoeff114Expr (.pow (.var 8) 4))) (.var 6)) (.mul ThetaCoeff015Expr (.pow (.var 8) 5)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul ThetaCoeff600Expr (.var 6)) (.mul ThetaCoeff501Expr (.var 8))) (.var 6)) (.mul ThetaCoeff402Expr (.pow (.var 8) 2))) (.var 6)) (.mul ThetaCoeff303Expr (.pow (.var 8) 3))) (.var 6)) (.mul ThetaCoeff204Expr (.pow (.var 8) 4))) (.var 6)) (.mul ThetaCoeff105Expr (.pow (.var 8) 5))) (.var 6)) (.mul ThetaCoeff006Expr (.pow (.var 8) 6))))
def JE : Expr :=
  (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul JCoeff090Expr (.var 7)) (.add (.mul JCoeff180Expr (.var 6)) (.mul JCoeff081Expr (.var 8)))) (.var 7)) (.add (.mul (.add (.mul JCoeff270Expr (.var 6)) (.mul JCoeff171Expr (.var 8))) (.var 6)) (.mul JCoeff072Expr (.pow (.var 8) 2)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul JCoeff360Expr (.var 6)) (.mul JCoeff261Expr (.var 8))) (.var 6)) (.mul JCoeff162Expr (.pow (.var 8) 2))) (.var 6)) (.mul JCoeff063Expr (.pow (.var 8) 3)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul (.add (.mul JCoeff450Expr (.var 6)) (.mul JCoeff351Expr (.var 8))) (.var 6)) (.mul JCoeff252Expr (.pow (.var 8) 2))) (.var 6)) (.mul JCoeff153Expr (.pow (.var 8) 3))) (.var 6)) (.mul JCoeff054Expr (.pow (.var 8) 4)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul JCoeff540Expr (.var 6)) (.mul JCoeff441Expr (.var 8))) (.var 6)) (.mul JCoeff342Expr (.pow (.var 8) 2))) (.var 6)) (.mul JCoeff243Expr (.pow (.var 8) 3))) (.var 6)) (.mul JCoeff144Expr (.pow (.var 8) 4))) (.var 6)) (.mul JCoeff045Expr (.pow (.var 8) 5)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul JCoeff630Expr (.var 6)) (.mul JCoeff531Expr (.var 8))) (.var 6)) (.mul JCoeff432Expr (.pow (.var 8) 2))) (.var 6)) (.mul JCoeff333Expr (.pow (.var 8) 3))) (.var 6)) (.mul JCoeff234Expr (.pow (.var 8) 4))) (.var 6)) (.mul JCoeff135Expr (.pow (.var 8) 5))) (.var 6)) (.mul JCoeff036Expr (.pow (.var 8) 6)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul JCoeff720Expr (.var 6)) (.mul JCoeff621Expr (.var 8))) (.var 6)) (.mul JCoeff522Expr (.pow (.var 8) 2))) (.var 6)) (.mul JCoeff423Expr (.pow (.var 8) 3))) (.var 6)) (.mul JCoeff324Expr (.pow (.var 8) 4))) (.var 6)) (.mul JCoeff225Expr (.pow (.var 8) 5))) (.var 6)) (.mul JCoeff126Expr (.pow (.var 8) 6))) (.var 6)) (.mul JCoeff027Expr (.pow (.var 8) 7)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul JCoeff810Expr (.var 6)) (.mul JCoeff711Expr (.var 8))) (.var 6)) (.mul JCoeff612Expr (.pow (.var 8) 2))) (.var 6)) (.mul JCoeff513Expr (.pow (.var 8) 3))) (.var 6)) (.mul JCoeff414Expr (.pow (.var 8) 4))) (.var 6)) (.mul JCoeff315Expr (.pow (.var 8) 5))) (.var 6)) (.mul JCoeff216Expr (.pow (.var 8) 6))) (.var 6)) (.mul JCoeff117Expr (.pow (.var 8) 7))) (.var 6)) (.mul JCoeff018Expr (.pow (.var 8) 8)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul JCoeff900Expr (.var 6)) (.mul JCoeff801Expr (.var 8))) (.var 6)) (.mul JCoeff702Expr (.pow (.var 8) 2))) (.var 6)) (.mul JCoeff603Expr (.pow (.var 8) 3))) (.var 6)) (.mul JCoeff504Expr (.pow (.var 8) 4))) (.var 6)) (.mul JCoeff405Expr (.pow (.var 8) 5))) (.var 6)) (.mul JCoeff306Expr (.pow (.var 8) 6))) (.var 6)) (.mul JCoeff207Expr (.pow (.var 8) 7))) (.var 6)) (.mul JCoeff108Expr (.pow (.var 8) 8))) (.var 6)) (.mul JCoeff009Expr (.pow (.var 8) 9))))

def h2 : Expr := .mul HE HE
def h3 : Expr := .mul h2 HE
def h4 : Expr := .mul h3 HE
def h5 : Expr := .mul h4 HE
def h6 : Expr := .mul h5 HE
def u2 : Expr := .mul UE UE
def u3 : Expr := .mul u2 UE
def u4 : Expr := .mul u3 UE
def q2 : Expr := .mul ThE ThE
def q3 : Expr := .mul q2 ThE
def s2 : Expr := .mul SExpr SExpr
def s3 : Expr := .mul s2 SExpr
def s4 : Expr := .mul s3 SExpr
def t2 : Expr := .mul TExpr TExpr
def LHS : Expr := .mul JE JE
def RHS : Expr :=
  (.add (.sub (.add (.mul (.num 4) q3) (.mul (.mul (.mul (.num 4) SExpr) ThE) h4)) (.mul TExpr h6))
  (.mul UE (.sub (.add (.sub (.sub (.sub (.sub (.sub (.add (.sub (.mul (.mul (.mul (.num 27) TExpr) UE) q2) (.mul (.mul (.mul (.num 108) s3) ThE) u3)) (.mul (.mul (.mul (.mul (.mul (.num 54) SExpr) TExpr) ThE) u2) HE)) (.mul (.mul (.mul (.mul (.num 72) s2) ThE) UE) h2)) (.mul (.mul (.mul (.num 18) TExpr) ThE) h3)) (.mul (.mul (.mul (.num 432) s4) u4) HE)) (.mul (.mul (.mul (.mul (.num 297) s2) TExpr) u3) h2)) (.mul (.mul (.mul (.num 108) t2) u2) h3)) (.mul (.mul (.mul (.mul (.num 54) SExpr) TExpr) UE) h4)) (.mul (.mul (.num 16) s2) h5))))

def base : ℕ := 2^50

theorem check : certify base LHS RHS=true := by decide +kernel

variable {K : Type*} [CommRing K]

section Evaluation
variable (a b c d e f x y z : K)
lemma vars : Lean.RArray.get (context a b c d e f x y z) 0=a ∧ Lean.RArray.get (context a b c d e f x y z) 1=b ∧ Lean.RArray.get (context a b c d e f x y z) 2=c ∧ Lean.RArray.get (context a b c d e f x y z) 3=d ∧ Lean.RArray.get (context a b c d e f x y z) 4=e ∧ Lean.RArray.get (context a b c d e f x y z) 5=f ∧ Lean.RArray.get (context a b c d e f x y z) 6=x ∧ Lean.RArray.get (context a b c d e f x y z) 7=y ∧ Lean.RArray.get (context a b c d e f x y z) 8=z :=
  ⟨rfl,rfl,rfl,rfl,rfl,rfl,rfl,rfl,rfl⟩
lemma U_eval : ringEval (fun i => (context a b c d e f x y z).get i) UE=CycleCubic.cubic a b c d e f x y z := by
  obtain ⟨h0,h1,h2,h3,h4,h5,h6,h7,h8⟩ := vars a b c d e f x y z
  simp only [UE,ringEval,h0,h1,h2,h3,h4,h5,h6,h7,h8,CycleCubic.cubic]; ring
lemma H_eval : ringEval (fun i => (context a b c d e f x y z).get i) HE=H a b c d e f x y z := by
  obtain ⟨-,-,-,-,-,-,h6,h7,h8⟩ := vars a b c d e f x y z
  simp only [HE,ringEval,h6,h7,h8]
  simp only [ringEval_denote,HCoeff300_denote,HCoeff210_denote,HCoeff201_denote,HCoeff120_denote,HCoeff111_denote,HCoeff102_denote,HCoeff030_denote,HCoeff021_denote,HCoeff012_denote,HCoeff003_denote,H]; ring
lemma Theta_eval : ringEval (fun i => (context a b c d e f x y z).get i) ThE=Theta a b c d e f x y z := by
  obtain ⟨-,-,-,-,-,-,h6,h7,h8⟩ := vars a b c d e f x y z
  simp only [ThE,ringEval,h6,h7,h8]
  simp only [ringEval_denote,ThetaCoeff600_denote,ThetaCoeff510_denote,ThetaCoeff501_denote,ThetaCoeff420_denote,ThetaCoeff411_denote,ThetaCoeff402_denote,ThetaCoeff330_denote,ThetaCoeff321_denote,ThetaCoeff312_denote,ThetaCoeff303_denote,ThetaCoeff240_denote,ThetaCoeff231_denote,ThetaCoeff222_denote,ThetaCoeff213_denote,ThetaCoeff204_denote,ThetaCoeff150_denote,ThetaCoeff141_denote,ThetaCoeff132_denote,ThetaCoeff123_denote,ThetaCoeff114_denote,ThetaCoeff105_denote,ThetaCoeff060_denote,ThetaCoeff051_denote,ThetaCoeff042_denote,ThetaCoeff033_denote,ThetaCoeff024_denote,ThetaCoeff015_denote,ThetaCoeff006_denote,Theta]; ring
lemma J_eval : ringEval (fun i => (context a b c d e f x y z).get i) JE=J a b c d e f x y z := by
  obtain ⟨-,-,-,-,-,-,h6,h7,h8⟩ := vars a b c d e f x y z
  simp only [JE,ringEval,h6,h7,h8]
  simp only [ringEval_denote,JCoeff900_denote,JCoeff810_denote,JCoeff801_denote,JCoeff720_denote,JCoeff711_denote,JCoeff702_denote,JCoeff630_denote,JCoeff621_denote,JCoeff612_denote,JCoeff603_denote,JCoeff540_denote,JCoeff531_denote,JCoeff522_denote,JCoeff513_denote,JCoeff504_denote,JCoeff450_denote,JCoeff441_denote,JCoeff432_denote,JCoeff423_denote,JCoeff414_denote,JCoeff405_denote,JCoeff360_denote,JCoeff351_denote,JCoeff342_denote,JCoeff333_denote,JCoeff324_denote,JCoeff315_denote,JCoeff306_denote,JCoeff270_denote,JCoeff261_denote,JCoeff252_denote,JCoeff243_denote,JCoeff234_denote,JCoeff225_denote,JCoeff216_denote,JCoeff207_denote,JCoeff180_denote,JCoeff171_denote,JCoeff162_denote,JCoeff153_denote,JCoeff144_denote,JCoeff135_denote,JCoeff126_denote,JCoeff117_denote,JCoeff108_denote,JCoeff090_denote,JCoeff081_denote,JCoeff072_denote,JCoeff063_denote,JCoeff054_denote,JCoeff045_denote,JCoeff036_denote,JCoeff027_denote,JCoeff018_denote,JCoeff009_denote,J]; ring
end Evaluation

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
  have h := certify_sound (fun i => ctx.get i) base LHS RHS check
  simp only [LHS,RHS,h2,h3,h4,h5,h6,u2,u3,u4,q2,q3,s2,s3,s4,t2,ringEval] at h
  simp only [ctx,U_eval,H_eval,Theta_eval,J_eval] at h
  simp only [ctx] at hs ht
  simp only [hs,ht,Int.cast_ofNat] at h
  linear_combination h
end CycleCovariants.Full
#print axioms CycleCovariants.Full.full_syzygy
