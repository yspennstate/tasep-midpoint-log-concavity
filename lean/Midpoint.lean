import Mathlib.Data.Int.Basic
import Mathlib.Tactic.Linarith

-- Integer balancing increases a two-factor product.
theorem balanced_product (a b r s : ℤ)
    (hsum : a + b = r + s)
    (hlo : -1 ≤ r - s) (hhi : r - s ≤ 1) :
    a * b ≤ r * s := by
  have hsq : (a + b)^2 = (r + s)^2 := congrArg (fun x : ℤ => x^2) hsum
  have hrs : r - s = -1 ∨ r - s = 0 ∨ r - s = 1 := by omega
  rcases lt_trichotomy a b with hab | hab | hab
  · have hgap : a + 1 ≤ b := by omega
    have hnonneg : 0 ≤ (b - a - 1)^2 := sq_nonneg _
    rcases hrs with hrs | hrs | hrs <;> nlinarith
  · have heq : r = s := by omega
    nlinarith
  · have hgap : b + 1 ≤ a := by omega
    have hnonneg : 0 ≤ (a - b - 1)^2 := sq_nonneg _
    rcases hrs with hrs | hrs | hrs <;> nlinarith

-- Difference factors improve under coordinate rounding, including integer shifts.
theorem rounded_difference_product (x₁ x₂ y₁ y₂ c : ℤ) :
    (x₁ - x₂ + c) * (y₁ - y₂ + c) ≤
    ((x₁ + y₁) / 2 - (x₂ + y₂) / 2 + c) *
    ((x₁ + y₁ + 1) / 2 - (x₂ + y₂ + 1) / 2 + c) := by
  apply balanced_product <;> omega

-- Coordinate comparisons with an integer offset survive both roundings.
theorem rounded_constraint (a b c d k : ℤ)
    (h₁ : a + k ≤ b) (h₂ : c + k ≤ d) :
    (a + c) / 2 + k ≤ (b + d) / 2 ∧
    (a + c + 1) / 2 + k ≤ (b + d + 1) / 2 := by
  omega

#print axioms balanced_product
#print axioms rounded_difference_product
#print axioms rounded_constraint
