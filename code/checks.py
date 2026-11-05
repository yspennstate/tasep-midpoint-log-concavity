"""Checks for the rounded-midpoint log-concavity of step TASEP.

1. Exact rational check of the array inequality W(A)W(B) <= W(R)W(S), R and S the rounded
   midpoints, for every pair of arrays with top row bounded by a cap (n = 2, 3, 4), at t = 1
   (the common factor e^{-nt} t^{|lambda|} is the same on both sides).
2. The edge marginal of the array law against the TASEP law computed by uniformization of the
   generator (n = 2 at t = 2, n = 3 at t = 1/2), 60 significant digits.
3. The theorem itself: p(a)p(b) <= p(floor((a+b)/2)) p(ceil((a+b)/2)) for every pair of
   states in a box, n = 2 and n = 3, t = 1/2, 1, 2, with p from the array formula.
Needs mpmath. Output: checks_out.txt (stdout) and checks.json.
"""
import itertools, json, math, pathlib
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 60
out = {}


def lower_rows(top):
    """All Gelfand-Tsetlin patterns with the given top row (weakly decreasing)."""
    if len(top) == 1:
        yield (top,)
        return
    for lower in itertools.product(*(range(top[i + 1], top[i] + 1) for i in range(len(top) - 1))):
        for pat in lower_rows(lower):
            yield pat + (top,)


def top_weight(top):
    """prod_{i<j}(l_i - l_j) / prod l_i!, l_i = lambda_i + n - i."""
    n = len(top)
    l = [a + n - i - 1 for i, a in enumerate(top)]
    num = math.prod(l[i] - l[j] for i in range(n) for j in range(i + 1, n))
    return Fraction(num, math.prod(math.factorial(x) for x in l))


# 1. exact array inequality
rows = []
for n, cap in [(2, 6), (3, 3), (4, 2)]:
    W = {}
    for s in itertools.combinations_with_replacement(range(cap + 1), n):
        top = tuple(reversed(s))
        w = top_weight(top)
        for pat in lower_rows(top):
            W[tuple(x for r in pat for x in r)] = w
    count = 0
    for A, B in itertools.combinations_with_replacement(W, 2):
        R = tuple((x + y) // 2 for x, y in zip(A, B))
        S = tuple((x + y + 1) // 2 for x, y in zip(A, B))
        assert R in W and S in W
        assert W[A] * W[B] <= W[R] * W[S], (A, B)
        count += 1
    rows.append(dict(n=n, cap=cap, arrays=len(W), pairs=count))
    print(f"1. n={n} top cap {cap}: {len(W)} arrays, {count} pairs, inequality holds exactly")
out["array_inequality"] = rows


# edge marginal of the array law, by direct summation over arrays with top row <= T
def edge_law(n, t, box, T):
    t = mp.mpf(t)
    p = {}
    for s in itertools.combinations_with_replacement(range(T + 1), n):
        top = tuple(reversed(s))
        if top[n - 1] > box:
            continue
        w = mp.exp(-n * t) * t ** sum(top) * mp.mpf(top_weight(top).numerator) / top_weight(top).denominator
        if n == 2:
            a, b = top
            for u in range(b, min(a, box) + 1):
                p[(u, b)] = p.get((u, b), 0) + w
        elif n == 3:
            a, b, c = top
            for v in range(c, min(b, box) + 1):
                for u in range(v, min(a, box) + 1):
                    # number of lambda^2_1 in [max(u, b), a]
                    p[(u, v, c)] = p.get((u, v, c), 0) + w * (a - max(u, b) + 1)
    return p


# TASEP in the coordinates D_i = X_i + i, by uniformization (states D_1 >= ... >= D_n >= 0)
def tasep_law(n, t, box, K=200):
    t = mp.mpf(t)
    states = [tuple(reversed(s)) for s in itertools.combinations_with_replacement(range(box + 1), n)]
    idx = {s: i for i, s in enumerate(states)}
    moves = []
    for s in states:
        dest, rate = [], 0
        for i in range(n):
            if i == 0 or s[i] < s[i - 1]:
                rate += 1
                z = list(s)
                z[i] += 1
                if tuple(z) in idx:
                    dest.append(idx[tuple(z)])
        moves.append((dest, n - rate))
    v = [mp.mpf(0)] * len(states)
    v[idx[(0,) * n]] = mp.mpf(1)
    w = mp.exp(-n * t)
    p = [w * x for x in v]
    for k in range(1, K + 1):
        nxt = [mp.mpf(0)] * len(v)
        for j, x in enumerate(v):
            if x:
                dest, stay = moves[j]
                nxt[j] += x * stay / n
                for h in dest:
                    nxt[h] += x / n
        v = nxt
        w *= n * t / k
        for j, x in enumerate(v):
            p[j] += w * x
    return dict(zip(states, p))


# 2. representation check
rows = []
for n, t, box in [(2, 2, 10), (3, 0.5, 6)]:
    pe = edge_law(n, t, box, 60)
    pt = tasep_law(n, t, box)
    err = max(abs(pe[s] - pt[s]) / pt[s] for s in pt)
    rows.append(dict(n=n, t=t, states=len(pt), max_relative_difference=float(err)))
    print(f"2. n={n} t={t}: {len(pt)} states, array marginal vs generator, max relative difference {mp.nstr(err, 3)}")
out["representation"] = rows

# 3. the theorem on a box
rows = []
for n, box in [(2, 10), (3, 6)]:
    for t in (0.5, 1, 2):
        p = edge_law(n, t, box, 70)
        states = list(p)
        count, worst = 0, mp.inf
        for a, b in itertools.combinations_with_replacement(states, 2):
            r = tuple((x + y) // 2 for x, y in zip(a, b))
            s = tuple((x + y + 1) // 2 for x, y in zip(a, b))
            if {r, s} == {a, b}:
                continue
            ratio = p[r] * p[s] / (p[a] * p[b])
            assert ratio > 1 - mp.mpf(10) ** -30, (a, b, ratio)
            worst = min(worst, ratio)
            count += 1
        rows.append(dict(n=n, t=t, states=len(states), pairs=count, smallest_ratio=float(worst)))
        print(f"3. n={n} t={t}: {len(states)} states, {count} nontrivial pairs, smallest ratio {mp.nstr(worst, 8)}")
out["theorem_on_box"] = rows

pathlib.Path(__file__).with_name("checks.json").write_text(json.dumps(out, indent=1))
