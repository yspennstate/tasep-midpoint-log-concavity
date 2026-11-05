# Log-concavity of the step TASEP distribution at rounded midpoints

Yitzchak Shmalo

This repository contains the paper and the computations for the following result. Let X(t) = (X_1(t), ..., X_n(t)) be the positions of n particles of the continuous-time totally asymmetric simple exclusion process started from the step configuration X_i(0) = -i, and let p_t(x) = P(X(t) = x) for x in Z^n. Then for every n, every t >= 0 and all a, b in Z^n,

    p_t(a) p_t(b) <= p_t(floor((a+b)/2)) p_t(ceil((a+b)/2)),

with rounding taken coordinatewise. In particular p_t(a) p_t(b) <= p_t((a+b)/2)^2 whenever (a+b)/2 is a lattice point. This answers the question on log-concavity of the TASEP transition probability in Section 3 of the open problem list of the AIM workshop *All roads to KPZ universality* (L. Petrov, A. Saenz, E. Sorensen, H. T. Vu, March 2025). The same holds for the law of any subset of the particles, so each single particle position has a log-concave law.

The proof writes X(t) as the diagonal of the Borodin-Ferrari interlacing array, whose fixed-time law has the weight e^{-nt} t^{|lambda|} prod_{i<j}(l_i - l_j) / prod_i l_i!, l_i = lambda_i + n - i, depending only on the top row lambda. This weight satisfies the rounded-midpoint inequality factor by factor, and the discrete four-functions theorem of Halikias, Klartag and Slomka passes the inequality to coordinate marginals.

## Contents

- `paper/main.tex`, `paper/main.pdf`: the paper.
- `code/checks.py`: exact rational check of the array inequality for every pair of patterns below a cap (n = 2, 3, 4); the diagonal marginal of the array law against the TASEP law obtained by uniformization of the generator, at 60 digits; and a direct check of the theorem on boxes of states for n = 2, 3 and t = 1/2, 1, 2. Output in `code/checks_out.txt` and `code/checks.json`.
- `lean/Midpoint.lean`: the balanced-product inequality, its consequence for rounded differences, and the preservation of shifted coordinate comparisons under rounding, checked with Lean 4 and Mathlib.

## Running

    cd code
    pip install mpmath
    python checks.py > checks_out.txt

The Lean file compiles inside any Lean 4 project with Mathlib (`lake env lean Midpoint.lean`).

## License

Paper: CC BY 4.0. Code: MIT.
