import Proposal.CubicCovariants
import HodgeComputer.Constructions.HoldoutCover.PolynomialCalculus
import Proposal.FullSyzygy

/-! The universal relative differential identity on the affine cubic chart.
The integer covariants use the classical Fisher/SageMath normalization credited
in CubicCovariants. The relative identity is checked as one Kronecker integer identity (FullSyzygy).
Prepared with OpenAI Codex assistance for hodgepodgeconjecture. -/
noncomputable section
namespace CycleCovariants
open MvPolynomial
set_option maxRecDepth 4096
set_option maxHeartbeats 2000000
set_option backward.isDefEq.respectTransparency false
variable {K : Type*} [CommRing K]
abbrev Affine (K : Type*) [CommRing K] := MvPolynomial (Fin 2) K
local notation "x" => (X (R := K) (0:Fin 2))
local notation "y" => (X (R := K) (1:Fin 2))
def Hpoly (a b c d e f : K) : Affine K :=
  ((((C (HCoeff300 a b c d e f))*x^(3:ℕ) + (C (HCoeff210 a b c d e f))*x^(2:ℕ)*y) + ((C (HCoeff201 a b c d e f))*x^(2:ℕ) + ((C (HCoeff120 a b c d e f))*x*y^(2:ℕ) + (C (HCoeff111 a b c d e f))*x*y))) + (((C (HCoeff102 a b c d e f))*x + (C (HCoeff030 a b c d e f))*y^(3:ℕ)) + ((C (HCoeff021 a b c d e f))*y^(2:ℕ) + ((C (HCoeff012 a b c d e f))*y + (C (HCoeff003 a b c d e f))*1))))
def Thetapoly (a b c d e f : K) : Affine K :=
  (((((C (ThetaCoeff600 a b c d e f))*x^(6:ℕ) + ((C (ThetaCoeff510 a b c d e f))*x^(5:ℕ)*y + (C (ThetaCoeff501 a b c d e f))*x^(5:ℕ))) + (((C (ThetaCoeff420 a b c d e f))*x^(4:ℕ)*y^(2:ℕ) + (C (ThetaCoeff411 a b c d e f))*x^(4:ℕ)*y) + ((C (ThetaCoeff402 a b c d e f))*x^(4:ℕ) + (C (ThetaCoeff330 a b c d e f))*x^(3:ℕ)*y^(3:ℕ)))) + (((C (ThetaCoeff321 a b c d e f))*x^(3:ℕ)*y^(2:ℕ) + ((C (ThetaCoeff312 a b c d e f))*x^(3:ℕ)*y + (C (ThetaCoeff303 a b c d e f))*x^(3:ℕ))) + (((C (ThetaCoeff240 a b c d e f))*x^(2:ℕ)*y^(4:ℕ) + (C (ThetaCoeff231 a b c d e f))*x^(2:ℕ)*y^(3:ℕ)) + ((C (ThetaCoeff222 a b c d e f))*x^(2:ℕ)*y^(2:ℕ) + (C (ThetaCoeff213 a b c d e f))*x^(2:ℕ)*y)))) + ((((C (ThetaCoeff204 a b c d e f))*x^(2:ℕ) + ((C (ThetaCoeff150 a b c d e f))*x*y^(5:ℕ) + (C (ThetaCoeff141 a b c d e f))*x*y^(4:ℕ))) + (((C (ThetaCoeff132 a b c d e f))*x*y^(3:ℕ) + (C (ThetaCoeff123 a b c d e f))*x*y^(2:ℕ)) + ((C (ThetaCoeff114 a b c d e f))*x*y + (C (ThetaCoeff105 a b c d e f))*x))) + (((C (ThetaCoeff060 a b c d e f))*y^(6:ℕ) + ((C (ThetaCoeff051 a b c d e f))*y^(5:ℕ) + (C (ThetaCoeff042 a b c d e f))*y^(4:ℕ))) + (((C (ThetaCoeff033 a b c d e f))*y^(3:ℕ) + (C (ThetaCoeff024 a b c d e f))*y^(2:ℕ)) + ((C (ThetaCoeff015 a b c d e f))*y + (C (ThetaCoeff006 a b c d e f))*1)))))
def Jpoly (a b c d e f : K) : Affine K :=
  ((((((C (JCoeff900 a b c d e f))*x^(9:ℕ) + ((C (JCoeff810 a b c d e f))*x^(8:ℕ)*y + (C (JCoeff801 a b c d e f))*x^(8:ℕ))) + ((C (JCoeff720 a b c d e f))*x^(7:ℕ)*y^(2:ℕ) + ((C (JCoeff711 a b c d e f))*x^(7:ℕ)*y + (C (JCoeff702 a b c d e f))*x^(7:ℕ)))) + (((C (JCoeff630 a b c d e f))*x^(6:ℕ)*y^(3:ℕ) + ((C (JCoeff621 a b c d e f))*x^(6:ℕ)*y^(2:ℕ) + (C (JCoeff612 a b c d e f))*x^(6:ℕ)*y)) + (((C (JCoeff603 a b c d e f))*x^(6:ℕ) + (C (JCoeff540 a b c d e f))*x^(5:ℕ)*y^(4:ℕ)) + ((C (JCoeff531 a b c d e f))*x^(5:ℕ)*y^(3:ℕ) + (C (JCoeff522 a b c d e f))*x^(5:ℕ)*y^(2:ℕ))))) + ((((C (JCoeff513 a b c d e f))*x^(5:ℕ)*y + ((C (JCoeff504 a b c d e f))*x^(5:ℕ) + (C (JCoeff450 a b c d e f))*x^(4:ℕ)*y^(5:ℕ))) + (((C (JCoeff441 a b c d e f))*x^(4:ℕ)*y^(4:ℕ) + (C (JCoeff432 a b c d e f))*x^(4:ℕ)*y^(3:ℕ)) + ((C (JCoeff423 a b c d e f))*x^(4:ℕ)*y^(2:ℕ) + (C (JCoeff414 a b c d e f))*x^(4:ℕ)*y))) + (((C (JCoeff405 a b c d e f))*x^(4:ℕ) + ((C (JCoeff360 a b c d e f))*x^(3:ℕ)*y^(6:ℕ) + (C (JCoeff351 a b c d e f))*x^(3:ℕ)*y^(5:ℕ))) + (((C (JCoeff342 a b c d e f))*x^(3:ℕ)*y^(4:ℕ) + (C (JCoeff333 a b c d e f))*x^(3:ℕ)*y^(3:ℕ)) + ((C (JCoeff324 a b c d e f))*x^(3:ℕ)*y^(2:ℕ) + (C (JCoeff315 a b c d e f))*x^(3:ℕ)*y))))) + (((((C (JCoeff306 a b c d e f))*x^(3:ℕ) + ((C (JCoeff270 a b c d e f))*x^(2:ℕ)*y^(7:ℕ) + (C (JCoeff261 a b c d e f))*x^(2:ℕ)*y^(6:ℕ))) + (((C (JCoeff252 a b c d e f))*x^(2:ℕ)*y^(5:ℕ) + (C (JCoeff243 a b c d e f))*x^(2:ℕ)*y^(4:ℕ)) + ((C (JCoeff234 a b c d e f))*x^(2:ℕ)*y^(3:ℕ) + (C (JCoeff225 a b c d e f))*x^(2:ℕ)*y^(2:ℕ)))) + (((C (JCoeff216 a b c d e f))*x^(2:ℕ)*y + ((C (JCoeff207 a b c d e f))*x^(2:ℕ) + (C (JCoeff180 a b c d e f))*x*y^(8:ℕ))) + (((C (JCoeff171 a b c d e f))*x*y^(7:ℕ) + (C (JCoeff162 a b c d e f))*x*y^(6:ℕ)) + ((C (JCoeff153 a b c d e f))*x*y^(5:ℕ) + (C (JCoeff144 a b c d e f))*x*y^(4:ℕ))))) + ((((C (JCoeff135 a b c d e f))*x*y^(3:ℕ) + ((C (JCoeff126 a b c d e f))*x*y^(2:ℕ) + (C (JCoeff117 a b c d e f))*x*y)) + (((C (JCoeff108 a b c d e f))*x + (C (JCoeff090 a b c d e f))*y^(9:ℕ)) + ((C (JCoeff081 a b c d e f))*y^(8:ℕ) + (C (JCoeff072 a b c d e f))*y^(7:ℕ)))) + (((C (JCoeff063 a b c d e f))*y^(6:ℕ) + ((C (JCoeff054 a b c d e f))*y^(5:ℕ) + (C (JCoeff045 a b c d e f))*y^(4:ℕ))) + (((C (JCoeff036 a b c d e f))*y^(3:ℕ) + (C (JCoeff027 a b c d e f))*y^(2:ℕ)) + ((C (JCoeff018 a b c d e f))*y + (C (JCoeff009 a b c d e f))*1))))))
def Upoly (a b c d e f : K) : Affine K :=
  C a*x^3+C b*y^3+C c+C d*x^2*y+C e*y^2+C f*x

section Certificate
open Lean.Grind.CommRing CycleSignedExpressions CycleNineKronecker CycleCovariants.Reflection
def UxE : Expr :=
  (.add (.mul (.mul (.mul (.num 2) (.var 3)) (.var 6)) (.var 7)) (.add (.mul (.mul (.mul (.num 3) (.var 0)) (.var 6)) (.var 6)) (.mul (.var 5) (.pow (.var 8) 2))))
def UyE : Expr :=
  (.add (.mul (.add (.mul (.mul (.num 3) (.var 1)) (.var 7)) (.mul (.mul (.num 2) (.var 4)) (.var 8))) (.var 7)) (.mul (.mul (.var 3) (.var 6)) (.var 6)))
def HxE : Expr :=
  (.add (.mul (.add (.mul HCoeff120Expr (.var 7)) (.add (.mul (.mul (.num 2) HCoeff210Expr) (.var 6)) (.mul HCoeff111Expr (.var 8)))) (.var 7)) (.add (.mul (.add (.mul (.mul (.num 3) HCoeff300Expr) (.var 6)) (.mul (.mul (.num 2) HCoeff201Expr) (.var 8))) (.var 6)) (.mul HCoeff102Expr (.pow (.var 8) 2))))
def HyE : Expr :=
  (.add (.mul (.add (.mul (.mul (.num 3) HCoeff030Expr) (.var 7)) (.add (.mul (.mul (.num 2) HCoeff120Expr) (.var 6)) (.mul (.mul (.num 2) HCoeff021Expr) (.var 8)))) (.var 7)) (.add (.mul (.add (.mul HCoeff210Expr (.var 6)) (.mul HCoeff111Expr (.var 8))) (.var 6)) (.mul HCoeff012Expr (.pow (.var 8) 2))))
def ThxE : Expr :=
  (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul ThetaCoeff150Expr (.var 7)) (.add (.mul (.mul (.num 2) ThetaCoeff240Expr) (.var 6)) (.mul ThetaCoeff141Expr (.var 8)))) (.var 7)) (.add (.mul (.add (.mul (.mul (.num 3) ThetaCoeff330Expr) (.var 6)) (.mul (.mul (.num 2) ThetaCoeff231Expr) (.var 8))) (.var 6)) (.mul ThetaCoeff132Expr (.pow (.var 8) 2)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul (.mul (.num 4) ThetaCoeff420Expr) (.var 6)) (.mul (.mul (.num 3) ThetaCoeff321Expr) (.var 8))) (.var 6)) (.mul (.mul (.num 2) ThetaCoeff222Expr) (.pow (.var 8) 2))) (.var 6)) (.mul ThetaCoeff123Expr (.pow (.var 8) 3)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.mul (.num 5) ThetaCoeff510Expr) (.var 6)) (.mul (.mul (.num 4) ThetaCoeff411Expr) (.var 8))) (.var 6)) (.mul (.mul (.num 3) ThetaCoeff312Expr) (.pow (.var 8) 2))) (.var 6)) (.mul (.mul (.num 2) ThetaCoeff213Expr) (.pow (.var 8) 3))) (.var 6)) (.mul ThetaCoeff114Expr (.pow (.var 8) 4)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.mul (.num 6) ThetaCoeff600Expr) (.var 6)) (.mul (.mul (.num 5) ThetaCoeff501Expr) (.var 8))) (.var 6)) (.mul (.mul (.num 4) ThetaCoeff402Expr) (.pow (.var 8) 2))) (.var 6)) (.mul (.mul (.num 3) ThetaCoeff303Expr) (.pow (.var 8) 3))) (.var 6)) (.mul (.mul (.num 2) ThetaCoeff204Expr) (.pow (.var 8) 4))) (.var 6)) (.mul ThetaCoeff105Expr (.pow (.var 8) 5))))
def ThyE : Expr :=
  (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.mul (.num 6) ThetaCoeff060Expr) (.var 7)) (.add (.mul (.mul (.num 5) ThetaCoeff150Expr) (.var 6)) (.mul (.mul (.num 5) ThetaCoeff051Expr) (.var 8)))) (.var 7)) (.add (.mul (.add (.mul (.mul (.num 4) ThetaCoeff240Expr) (.var 6)) (.mul (.mul (.num 4) ThetaCoeff141Expr) (.var 8))) (.var 6)) (.mul (.mul (.num 4) ThetaCoeff042Expr) (.pow (.var 8) 2)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul (.mul (.num 3) ThetaCoeff330Expr) (.var 6)) (.mul (.mul (.num 3) ThetaCoeff231Expr) (.var 8))) (.var 6)) (.mul (.mul (.num 3) ThetaCoeff132Expr) (.pow (.var 8) 2))) (.var 6)) (.mul (.mul (.num 3) ThetaCoeff033Expr) (.pow (.var 8) 3)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.mul (.num 2) ThetaCoeff420Expr) (.var 6)) (.mul (.mul (.num 2) ThetaCoeff321Expr) (.var 8))) (.var 6)) (.mul (.mul (.num 2) ThetaCoeff222Expr) (.pow (.var 8) 2))) (.var 6)) (.mul (.mul (.num 2) ThetaCoeff123Expr) (.pow (.var 8) 3))) (.var 6)) (.mul (.mul (.num 2) ThetaCoeff024Expr) (.pow (.var 8) 4)))) (.var 7)) (.add (.mul (.add (.mul (.add (.mul (.add (.mul (.add (.mul ThetaCoeff510Expr (.var 6)) (.mul ThetaCoeff411Expr (.var 8))) (.var 6)) (.mul ThetaCoeff312Expr (.pow (.var 8) 2))) (.var 6)) (.mul ThetaCoeff213Expr (.pow (.var 8) 3))) (.var 6)) (.mul ThetaCoeff114Expr (.pow (.var 8) 4))) (.var 6)) (.mul ThetaCoeff015Expr (.pow (.var 8) 5))))

/-- Homogeneous form of the identity: the Jacobian term carries the factor `z`. -/
def DL : Expr :=
  .add (.sub (.mul Full.HE (.sub (.mul UxE ThyE) (.mul UyE ThxE)))
      (.mul (.mul (.num 2) Full.ThE) (.sub (.mul UxE HyE) (.mul UyE HxE))))
    (.mul (.mul (.num 3) (.var 8)) Full.JE)
def DR : Expr := .mul Full.UE (.sub (.mul HxE ThyE) (.mul HyE ThxE))

theorem dcheck : certify Full.base DL DR=true := by decide +kernel

variable (a b c d e f : K)
lemma gets : Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) 0=C a ∧ Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) 1=C b ∧
    Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) 2=C c ∧ Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) 3=C d ∧
    Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) 4=C e ∧ Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) 5=C f ∧
    Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) 6=x ∧ Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) 7=y ∧
    Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) 8=1 := ⟨rfl,rfl,rfl,rfl,rfl,rfl,rfl,rfl,rfl⟩

lemma U_w : ringEval (fun i => Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) i) Full.UE=Upoly a b c d e f := by
  rw [Full.U_eval]; simp only [CycleCubic.cubic,Upoly]; ring
lemma H_w : ringEval (fun i => Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) i) Full.HE=Hpoly a b c d e f := by
  rw [Full.H_eval]
  simp only [H,Hpoly,HCoeff300,HCoeff210,HCoeff201,HCoeff120,HCoeff111,HCoeff102,HCoeff030,HCoeff021,HCoeff012,HCoeff003,map_add,map_sub,map_mul,map_pow,map_neg,map_ofNat]; ring
lemma Th_w : ringEval (fun i => Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) i) Full.ThE=Thetapoly a b c d e f := by
  rw [Full.Theta_eval]
  simp only [Theta,Thetapoly,ThetaCoeff600,ThetaCoeff510,ThetaCoeff501,ThetaCoeff420,ThetaCoeff411,ThetaCoeff402,ThetaCoeff330,ThetaCoeff321,ThetaCoeff312,ThetaCoeff303,ThetaCoeff240,ThetaCoeff231,ThetaCoeff222,ThetaCoeff213,ThetaCoeff204,ThetaCoeff150,ThetaCoeff141,ThetaCoeff132,ThetaCoeff123,ThetaCoeff114,ThetaCoeff105,ThetaCoeff060,ThetaCoeff051,ThetaCoeff042,ThetaCoeff033,ThetaCoeff024,ThetaCoeff015,ThetaCoeff006,map_add,map_sub,map_mul,map_pow,map_neg,map_ofNat]; ring
lemma J_w : ringEval (fun i => Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) i) Full.JE=Jpoly a b c d e f := by
  rw [Full.J_eval]
  simp only [J,Jpoly,JCoeff900,JCoeff810,JCoeff801,JCoeff720,JCoeff711,JCoeff702,JCoeff630,JCoeff621,JCoeff612,JCoeff603,JCoeff540,JCoeff531,JCoeff522,JCoeff513,JCoeff504,JCoeff450,JCoeff441,JCoeff432,JCoeff423,JCoeff414,JCoeff405,JCoeff360,JCoeff351,JCoeff342,JCoeff333,JCoeff324,JCoeff315,JCoeff306,JCoeff270,JCoeff261,JCoeff252,JCoeff243,JCoeff234,JCoeff225,JCoeff216,JCoeff207,JCoeff180,JCoeff171,JCoeff162,JCoeff153,JCoeff144,JCoeff135,JCoeff126,JCoeff117,JCoeff108,JCoeff090,JCoeff081,JCoeff072,JCoeff063,JCoeff054,JCoeff045,JCoeff036,JCoeff027,JCoeff018,JCoeff009,map_add,map_sub,map_mul,map_pow,map_neg,map_ofNat]; ring
lemma Ux_w : ringEval (fun i => Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) i) UxE=pderiv 0 (Upoly a b c d e f) := by
  obtain ⟨g0,g1,g2,g3,g4,g5,g6,g7,g8⟩ := gets a b c d e f
  simp only [UxE,ringEval,g0,g1,g2,g3,g4,g5,g6,g7,g8]
  simp only [Upoly,map_add,map_neg,map_sub,pderiv_mul,pderiv_pow,pderiv_C,pderiv_X,Pi.single_apply]
  norm_num
  ring
lemma Uy_w : ringEval (fun i => Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) i) UyE=pderiv 1 (Upoly a b c d e f) := by
  obtain ⟨g0,g1,g2,g3,g4,g5,g6,g7,g8⟩ := gets a b c d e f
  simp only [UyE,ringEval,g0,g1,g2,g3,g4,g5,g6,g7,g8]
  simp only [Upoly,map_add,map_neg,map_sub,pderiv_mul,pderiv_pow,pderiv_C,pderiv_X,Pi.single_apply]
  norm_num
  ring
lemma Hx_w : ringEval (fun i => Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) i) HxE=pderiv 0 (Hpoly a b c d e f) := by
  obtain ⟨g0,g1,g2,g3,g4,g5,g6,g7,g8⟩ := gets a b c d e f
  simp only [HxE,ringEval,g0,g1,g2,g3,g4,g5,g6,g7,g8]
  simp only [ringEval_denote,HCoeff300_denote,HCoeff210_denote,HCoeff201_denote,HCoeff120_denote,HCoeff111_denote,HCoeff102_denote,HCoeff030_denote,HCoeff021_denote,HCoeff012_denote,HCoeff003_denote,Hpoly,map_add,map_neg,map_sub,pderiv_mul,pderiv_pow,pderiv_C,pderiv_X,Pi.single_apply]
  norm_num
  simp only [HCoeff300,HCoeff210,HCoeff201,HCoeff120,HCoeff111,HCoeff102,HCoeff030,HCoeff021,HCoeff012,HCoeff003,map_add,map_sub,map_mul,map_pow,map_neg,map_ofNat]
  ring
lemma Hy_w : ringEval (fun i => Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) i) HyE=pderiv 1 (Hpoly a b c d e f) := by
  obtain ⟨g0,g1,g2,g3,g4,g5,g6,g7,g8⟩ := gets a b c d e f
  simp only [HyE,ringEval,g0,g1,g2,g3,g4,g5,g6,g7,g8]
  simp only [ringEval_denote,HCoeff300_denote,HCoeff210_denote,HCoeff201_denote,HCoeff120_denote,HCoeff111_denote,HCoeff102_denote,HCoeff030_denote,HCoeff021_denote,HCoeff012_denote,HCoeff003_denote,Hpoly,map_add,map_neg,map_sub,pderiv_mul,pderiv_pow,pderiv_C,pderiv_X,Pi.single_apply]
  norm_num
  simp only [HCoeff300,HCoeff210,HCoeff201,HCoeff120,HCoeff111,HCoeff102,HCoeff030,HCoeff021,HCoeff012,HCoeff003,map_add,map_sub,map_mul,map_pow,map_neg,map_ofNat]
  ring
lemma Thx_w : ringEval (fun i => Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) i) ThxE=pderiv 0 (Thetapoly a b c d e f) := by
  obtain ⟨g0,g1,g2,g3,g4,g5,g6,g7,g8⟩ := gets a b c d e f
  simp only [ThxE,ringEval,g0,g1,g2,g3,g4,g5,g6,g7,g8]
  simp only [ringEval_denote,ThetaCoeff600_denote,ThetaCoeff510_denote,ThetaCoeff501_denote,ThetaCoeff420_denote,ThetaCoeff411_denote,ThetaCoeff402_denote,ThetaCoeff330_denote,ThetaCoeff321_denote,ThetaCoeff312_denote,ThetaCoeff303_denote,ThetaCoeff240_denote,ThetaCoeff231_denote,ThetaCoeff222_denote,ThetaCoeff213_denote,ThetaCoeff204_denote,ThetaCoeff150_denote,ThetaCoeff141_denote,ThetaCoeff132_denote,ThetaCoeff123_denote,ThetaCoeff114_denote,ThetaCoeff105_denote,ThetaCoeff060_denote,ThetaCoeff051_denote,ThetaCoeff042_denote,ThetaCoeff033_denote,ThetaCoeff024_denote,ThetaCoeff015_denote,ThetaCoeff006_denote,Thetapoly,map_add,map_neg,map_sub,pderiv_mul,pderiv_pow,pderiv_C,pderiv_X,Pi.single_apply]
  norm_num
  simp only [ThetaCoeff600,ThetaCoeff510,ThetaCoeff501,ThetaCoeff420,ThetaCoeff411,ThetaCoeff402,ThetaCoeff330,ThetaCoeff321,ThetaCoeff312,ThetaCoeff303,ThetaCoeff240,ThetaCoeff231,ThetaCoeff222,ThetaCoeff213,ThetaCoeff204,ThetaCoeff150,ThetaCoeff141,ThetaCoeff132,ThetaCoeff123,ThetaCoeff114,ThetaCoeff105,ThetaCoeff060,ThetaCoeff051,ThetaCoeff042,ThetaCoeff033,ThetaCoeff024,ThetaCoeff015,ThetaCoeff006,map_add,map_sub,map_mul,map_pow,map_neg,map_ofNat]
  ring
lemma Thy_w : ringEval (fun i => Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) i) ThyE=pderiv 1 (Thetapoly a b c d e f) := by
  obtain ⟨g0,g1,g2,g3,g4,g5,g6,g7,g8⟩ := gets a b c d e f
  simp only [ThyE,ringEval,g0,g1,g2,g3,g4,g5,g6,g7,g8]
  simp only [ringEval_denote,ThetaCoeff600_denote,ThetaCoeff510_denote,ThetaCoeff501_denote,ThetaCoeff420_denote,ThetaCoeff411_denote,ThetaCoeff402_denote,ThetaCoeff330_denote,ThetaCoeff321_denote,ThetaCoeff312_denote,ThetaCoeff303_denote,ThetaCoeff240_denote,ThetaCoeff231_denote,ThetaCoeff222_denote,ThetaCoeff213_denote,ThetaCoeff204_denote,ThetaCoeff150_denote,ThetaCoeff141_denote,ThetaCoeff132_denote,ThetaCoeff123_denote,ThetaCoeff114_denote,ThetaCoeff105_denote,ThetaCoeff060_denote,ThetaCoeff051_denote,ThetaCoeff042_denote,ThetaCoeff033_denote,ThetaCoeff024_denote,ThetaCoeff015_denote,ThetaCoeff006_denote,Thetapoly,map_add,map_neg,map_sub,pderiv_mul,pderiv_pow,pderiv_C,pderiv_X,Pi.single_apply]
  norm_num
  simp only [ThetaCoeff600,ThetaCoeff510,ThetaCoeff501,ThetaCoeff420,ThetaCoeff411,ThetaCoeff402,ThetaCoeff330,ThetaCoeff321,ThetaCoeff312,ThetaCoeff303,ThetaCoeff240,ThetaCoeff231,ThetaCoeff222,ThetaCoeff213,ThetaCoeff204,ThetaCoeff150,ThetaCoeff141,ThetaCoeff132,ThetaCoeff123,ThetaCoeff114,ThetaCoeff105,ThetaCoeff060,ThetaCoeff051,ThetaCoeff042,ThetaCoeff033,ThetaCoeff024,ThetaCoeff015,ThetaCoeff006,map_add,map_sub,map_mul,map_pow,map_neg,map_ofNat]
  ring
end Certificate

open CycleSignedExpressions CycleNineKronecker CycleCovariants.Reflection in
lemma relative_identity (a b c d e f : K) :
    let u := Upoly a b c d e f
    let h := Hpoly a b c d e f
    let q := Thetapoly a b c d e f
    let j := Jpoly a b c d e f
    h*(pderiv 0 u*pderiv 1 q-pderiv 1 u*pderiv 0 q)-
      2*q*(pderiv 0 u*pderiv 1 h-pderiv 1 u*pderiv 0 h)+3*j=
      u*(pderiv 0 h*pderiv 1 q-pderiv 1 h*pderiv 0 q) := by
  dsimp only
  have h := certify_sound (fun i => Lean.RArray.get (context (C a) (C b) (C c) (C d) (C e) (C f) x y (1 : Affine K)) i) Full.base DL DR dcheck
  obtain ⟨-,-,-,-,-,-,-,-,g8⟩ := gets a b c d e f
  simp only [DL,DR,ringEval,g8] at h
  rw [U_w,H_w,Th_w,J_w,Ux_w,Uy_w,Hx_w,Hy_w,Thx_w,Thy_w] at h
  simp only [Int.cast_ofNat,mul_one] at h
  linear_combination h
#print axioms relative_identity
end CycleCovariants
