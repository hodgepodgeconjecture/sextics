import Mathlib.RingTheory.Derivation.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Ring

/-!
Exact compact reduction for a cyclic cubic cover.

Research lemma for cyclic-triple-residual. This module does not assert
HasCertificate or assume existence of the missing scalar parameter current.
The determinant transport and incidence binding remain separate obligations.
Prepared with OpenAI Codex assistance for hodgepodgeconjecture.
-/

noncomputable section

namespace CyclicCubicReduction

variable {K : Type*} [Field K] [CharZero K]

def cubic (a b c d t : K) : K := a*t^3+b*t^2+c*t+d

def cubicDerivative (a b c t : K) : K := 3*a*t^2+2*b*t+c

def discriminant (a b c d : K) : K :=
  b^2*c^2-4*a*c^3-4*b^3*d-27*a^2*d^2+18*a*b*c*d

def discriminantJet (a b c d a₁ b₁ c₁ d₁ : K) : K :=
  (-54*a*d^2+18*b*c*d-4*c^3)*a₁ +
  (18*a*c*d-12*b^2*d+2*b*c^2)*b₁ +
  (18*a*b*d-12*a*c^2+2*b^2*c)*c₁ +
  (-54*a^2*d+18*a*b*c-4*b^3)*d₁

omit [CharZero K] in
/-- Bind the coefficient jet to differentiation of the actual cubic. -/
theorem cubic_derivation (E : Derivation ℤ K K) (a b c d t : K) :
    E (cubic a b c d t) = cubic (E a) (E b) (E c) (E d) t +
      cubicDerivative a b c t*E t := by
  simp [cubic, cubicDerivative, Derivation.leibniz, Derivation.leibniz_pow,
    smul_eq_mul]
  ring

omit [CharZero K] in
/-- Bind the discriminant jet to differentiation of the discriminant. -/
theorem discriminant_derivation (E : Derivation ℤ K K) (a b c d : K) :
    E (discriminant a b c d) = discriminantJet a b c d (E a) (E b) (E c) (E d) := by
  have h4 : E (4:K) = 0 := E.map_natCast 4
  have h18 : E (18:K) = 0 := E.map_natCast 18
  have h27 : E (27:K) = 0 := E.map_natCast 27
  simp [discriminant, discriminantJet, Derivation.leibniz,
    Derivation.leibniz_pow, smul_eq_mul, h4, h18, h27]
  ring

def numerator0 (a b c d a₁ b₁ c₁ d₁ : K) : K :=
  (6*b*d^2-2*c^2*d)*a₁ + (-9*a*d^2+b*c*d)*b₁ +
  (6*a*c*d-2*b^2*d)*c₁ + (3*a*b*d-4*a*c^2+b^2*c)*d₁

def numerator1 (a b c d a₁ b₁ c₁ d₁ : K) : K :=
  (-9*a*d^2+7*b*c*d-2*c^3)*a₁ + (-3*a*c*d-2*b^2*d+b*c^2)*b₁ +
  (3*a*b*d+2*a*c^2-b^2*c)*c₁ + (9*a^2*d-7*a*b*c+2*b^3)*d₁

def numerator2 (a b c d a₁ b₁ c₁ d₁ : K) : K :=
  (-3*a*c*d+4*b^2*d-b*c^2)*a₁ + (-6*a*b*d+2*a*c^2)*b₁ +
  (9*a^2*d-a*b*c)*c₁ + (-6*a^2*c+2*a*b^2)*d₁

omit [CharZero K] in
/-- The exact polynomial identity underlying the compact connection. -/
theorem compact_polynomial_identity (a b c d a₁ b₁ c₁ d₁ t : K) :
    let n₀ := numerator0 a b c d a₁ b₁ c₁ d₁
    let n₁ := numerator1 a b c d a₁ b₁ c₁ d₁
    let n₂ := numerator2 a b c d a₁ b₁ c₁ d₁
    2*(3*cubic a b c d t*(2*n₂*t+n₁) -
       2*cubicDerivative a b c t*(n₂*t^2+n₁*t+n₀)) =
      -4*discriminant a b c d*cubic a₁ b₁ c₁ d₁ t +
        discriminantJet a b c d a₁ b₁ c₁ d₁*cubic a b c d t := by
  dsimp [numerator0, numerator1, numerator2, cubic, cubicDerivative,
    discriminant, discriminantJet]
  ring

/-- Rational differentiation in the cover field; all denominators are explicit. -/
theorem cubic_primitive_rule (T : Derivation ℤ K K)
    (P P₁ w Z Z₁ : K) (hw : w ≠ 0) (hcover : w^3 = Z)
    (hP : T P = P₁) (hwjet : 3*w^2*T w = Z₁) :
    T (P/(3*w^2)) = (3*Z*P₁-2*Z₁*P)/(9*w^5) := by
  have hden : (3:K)*w^2 ≠ 0 := mul_ne_zero (by norm_num) (pow_ne_zero _ hw)
  have hDw : T w = Z₁/(3*w^2) := (eq_div_iff hden).2 (by simpa [mul_comm] using hwjet)
  have h3 : T (3:K) = 0 := T.map_natCast 3
  rw [Derivation.leibniz_div]
  simp only [Derivation.leibniz, Derivation.leibniz_pow, h3, hP, hDw,
    smul_eq_mul, nsmul_eq_mul]
  rw [← hcover]
  field_simp
  ring

/-- The residual scalar coefficient after removing the exact time derivative. -/
theorem compact_rational_identity
    (J J₁ w Z Z₁ Zs Δ Δs N N₁ : K)
    (hw : w ≠ 0) (hΔ : Δ ≠ 0) (hcover : w^3 = Z)
    (hpoly : 2*(3*Z*N₁-2*Z₁*N) = -4*Δ*Zs+Δs*Z) :
    (3*Z*J₁-2*J*Zs)/(9*w^5) =
      (3*Z*(J*N₁/Δ)-2*Z₁*(J*N/Δ))/(9*w^5) +
      (J₁/3-J*Δs/(18*Δ))/w^2 := by
  rw [← hcover] at hpoly ⊢
  field_simp
  linear_combination -27*J*hpoly

/-- Assemble the exact derivative and the remaining compact scalar. -/
theorem compact_source_of_quadratic (T : Derivation ℤ K K)
    (t J J₁ w Z Z₁ Zs Δ Δs n₀ n₁ n₂ : K)
    (hw : w ≠ 0) (hΔ : Δ ≠ 0) (hcover : w^3 = Z)
    (ht : T t = 1) (hJ : T J = 0) (hD : T Δ = 0)
    (h₀ : T n₀ = 0) (h₁ : T n₁ = 0) (h₂ : T n₂ = 0)
    (hwjet : 3*w^2*T w = Z₁)
    (hpoly : 2*(3*Z*(2*n₂*t+n₁)-2*Z₁*(n₂*t^2+n₁*t+n₀)) =
      -4*Δ*Zs+Δs*Z) :
    (3*Z*J₁-2*J*Zs)/(9*w^5) =
      T ((J*(n₂*t^2+n₁*t+n₀)/Δ)/(3*w^2)) +
      (J₁/3-J*Δs/(18*Δ))/w^2 := by
  have hP : T (J*(n₂*t^2+n₁*t+n₀)/Δ) = J*(2*n₂*t+n₁)/Δ := by
    simp only [Derivation.leibniz_div, Derivation.leibniz, Derivation.leibniz_pow,
      map_add, h₀, h₁, h₂, ht, hJ, hD, smul_eq_mul, nsmul_eq_mul]
    field_simp
    ring
  rw [cubic_primitive_rule T _ _ w Z Z₁ hw hcover hP hwjet]
  exact compact_rational_identity J J₁ w Z Z₁ Zs Δ Δs
    (n₂*t^2+n₁*t+n₀) (2*n₂*t+n₁) hw hΔ hcover hpoly

/-- The actual cubic identity supplies the polynomial premise of the reduction. -/
theorem compact_source (T : Derivation ℤ K K)
    (a b c d a₁ b₁ c₁ d₁ t J J₁ w : K)
    (hw : w ≠ 0) (hΔ : discriminant a b c d ≠ 0)
    (hcover : w^3 = cubic a b c d t)
    (ht : T t = 1) (hJ : T J = 0) (hD : T (discriminant a b c d) = 0)
    (h₀ : T (numerator0 a b c d a₁ b₁ c₁ d₁) = 0)
    (h₁ : T (numerator1 a b c d a₁ b₁ c₁ d₁) = 0)
    (h₂ : T (numerator2 a b c d a₁ b₁ c₁ d₁) = 0)
    (hwjet : 3*w^2*T w = cubicDerivative a b c t) :
    (3*cubic a b c d t*J₁-2*J*cubic a₁ b₁ c₁ d₁ t)/(9*w^5) =
      T ((J*(numerator2 a b c d a₁ b₁ c₁ d₁*t^2 +
        numerator1 a b c d a₁ b₁ c₁ d₁*t +
        numerator0 a b c d a₁ b₁ c₁ d₁)/discriminant a b c d)/(3*w^2)) +
      (J₁/3-J*discriminantJet a b c d a₁ b₁ c₁ d₁/
        (18*discriminant a b c d))/w^2 := by
  exact compact_source_of_quadratic T t J J₁ w
    (cubic a b c d t) (cubicDerivative a b c t) (cubic a₁ b₁ c₁ d₁ t)
    (discriminant a b c d) (discriminantJet a b c d a₁ b₁ c₁ d₁)
    (numerator0 a b c d a₁ b₁ c₁ d₁) (numerator1 a b c d a₁ b₁ c₁ d₁)
    (numerator2 a b c d a₁ b₁ c₁ d₁)
    hw hΔ hcover ht hJ hD h₀ h₁ h₂ hwjet
    (compact_polynomial_identity a b c d a₁ b₁ c₁ d₁ t)

/-- On the sixth-root cover the scalar connection is an ordinary derivative. -/
theorem kummer_derivative (E : Derivation ℤ K K) (b δ Δ : K)
    (hδ : δ ≠ 0) (hcover : δ^6 = Δ) :
    E (b/δ) = (E b - b*E Δ/(6*Δ))/δ := by
  rw [Derivation.leibniz_div, ← hcover]
  simp only [Derivation.leibniz_pow, smul_eq_mul, nsmul_eq_mul]
  field_simp [hδ]
  ring

open scoped BigOperators

omit [CharZero K] in
/-- A solution of the scalar equation reconstructs the parameter and time
currents. This theorem takes that equation as a premise; it does not solve it. -/
theorem current_reconstruction {ι : Type*} [Fintype ι]
    (T : Derivation ℤ K K) (E : ι → Derivation ℤ K K)
    (P B q : K) (b γ c S : ι → K)
    (hb : ∀ i, T (b i) = 0)
    (hconnection : ∀ i, E i q = γ i*q+T (c i))
    (hcurrent : ∑ i, (E i (b i)+γ i*b i) = B) :
    T P+B*q-∑ i, E i (S i) =
      T (P-∑ i, b i*c i)+∑ i, E i (b i*q-S i) := by
  have hterm (i : ι) : E i (b i*q-S i) =
      (E i (b i)+γ i*b i)*q+T (b i*c i)-E i (S i) := by
    simp only [map_sub, Derivation.leibniz, smul_eq_mul, hb, hconnection]
    ring
  simp only [hterm, map_sub, map_sum, Finset.sum_sub_distrib,
    Finset.sum_add_distrib, ← Finset.sum_mul, hcurrent]
  ring

end CyclicCubicReduction

#print axioms CyclicCubicReduction.compact_polynomial_identity
#print axioms CyclicCubicReduction.cubic_derivation
#print axioms CyclicCubicReduction.discriminant_derivation
#print axioms CyclicCubicReduction.cubic_primitive_rule
#print axioms CyclicCubicReduction.compact_rational_identity
#print axioms CyclicCubicReduction.compact_source_of_quadratic
#print axioms CyclicCubicReduction.compact_source
#print axioms CyclicCubicReduction.kummer_derivative
#print axioms CyclicCubicReduction.current_reconstruction
