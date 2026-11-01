"""Cross-check p_hunter_tf.tf_window against brute force and known large factors."""
import random
import sys
import time

import p_hunter_tf as tf

random.seed(1)
bad = 0


def brute(p, k_lo, k_hi):
    for k in range(max(k_lo, 1), k_hi + 1):
        q = 2 * k * p + 1
        if q >= (1 << p) - 1:
            return None
        if pow(2, p, q) == 1:
            return (q, k)
    return None


def check(label, got, want):
    global bad
    if got != want:
        bad += 1
        print(f"MISMATCH {label}: got={got} want={want}")


def small_primes(lo, hi):
    return [n for n in range(lo, hi) if n > 2 and all(n % d for d in range(2, int(n**0.5) + 1))]


# 1. small p, windows starting at 1 (includes q > 2^64 impossible here; covers min-window logic)
for p in small_primes(3, 400):
    check(f"p={p} full", tf.tf_window(p, 1, 20000), brute(p, 1, 20000))

# 2. random mid p, band-split must equal full scan
pool = small_primes(100_000, 101_500)
for p in random.sample(pool, 25):
    hi = 300_000
    full = brute(p, 1, hi)
    check(f"p={p} full", tf.tf_window(p, 1, hi), full)
    mid = hi // 3
    a = tf.tf_window(p, 1, mid)
    b = tf.tf_window(p, mid + 1, hi) if a is None else None
    check(f"p={p} split", a or b, full)

# 3. known factors above 2^64 (exercises the 128-bit path)
for p, q in ((163, 36230454570129675721),):
    k = (q - 1) // (2 * p)
    check(f"p={p} big q", tf.tf_window(p, k, k), (q, k))
    check(f"p={p} big q miss", tf.tf_window(p, k + 1, k + 1000), None)

# 4. known Mersenne primes must never be eliminated (window small enough to be fast)
for p in (521, 607, 1279, 2203, 4423, 9689, 21701, 44497, 86243, 132049,
          756839, 6972593, 13466917, 136279841):
    r = tf.tf_window(p, 1, 200_000)
    check(f"known prime p={p}", r, None)

# 5. big-p timing sanity
t0 = time.time()
r = tf.tf_window(136279843, 1, 50_000_000)
print(f"p=136279843 k<=5e7: {r} in {time.time()-t0:.2f}s")

print("OK" if not bad else f"{bad} FAILURES")
sys.exit(1 if bad else 0)
