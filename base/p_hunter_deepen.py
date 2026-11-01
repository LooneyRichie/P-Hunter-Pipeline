"""
p_hunter_deepen.py — iterative deepening B2 filter.

Round 1: B2 over k ∈ [1, K0] on the initial survivor set.
Round N: B2 over k ∈ (K_{N-2}, K_{N-1}] on round-(N-1) survivors only.

Each round the survivor set shrinks (fewer p's to probe) AND the k-work
per p is only the NEW band.  Total cost stays roughly
sum_i (survivors_i × k_step_i), not (survivors_0 × k_final).

Rust backend expects:
    p_hunter_rs.gate_b2_factor_probe(p, min_k, max_k) -> Optional[(q, k)]
Falls back to a pure-Python gate if the module isn't importable.
"""

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------
# Known constants
# ---------------------------------------------------------------

KNOWN_MERSENNE_EXPONENTS = {
    2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279,
    2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701,
    23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433,
    1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011,
    24036583, 25964951, 30402457, 32582657, 37156667, 42643801,
    43112609, 57885161, 74207281, 77232917, 82589933, 136279841,
}
KNOWN_COMPOSITE_DOUBLE_MERSENNE_EXPONENTS = {131071, 524287, 2147483647}
SMALL_KNOWN_P = frozenset({2, 3, 5, 7})
MR_BASES_64 = (2, 325, 9375, 28178, 450775, 9780504, 1795265022)

# ---------------------------------------------------------------
# Pure Python for base test testing against Rust backend
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Minimal filters (self-contained)
# ---------------------------------------------------------------

def is_prime_64(n):
    if n < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n == sp:  return True
        if n % sp == 0: return False
    d = n - 1; r = 0
    while d % 2 == 0:
        r += 1; d //= 2
    for base in MR_BASES_64:
        if base % n == 0: continue
        w = pow(base, d, n)
        if w == 1 or w == n - 1: continue
        for _ in range(r - 1):
            w = (w * w) % n
            if w == n - 1: break
        else:
            return False
    return True


_WHEEL_MOD = 2 * 3 * 5 * 7 * 11 * 13

def _build_wheel():
    w = bytearray(b"\x01") * _WHEEL_MOD
    for m in (2, 3, 5, 7, 11, 13):
        w[0::m] = b"\x00" * (_WHEEL_MOD // m)
    return bytes(w)

_WHEEL = _build_wheel()

def gate0(p):
    if p < 2: return False
    if p <= 13: return p in (2, 3, 5, 7, 11, 13)
    return _WHEEL[p % _WHEEL_MOD] == 1

# ---------------------------------------------------------------
# B2 dispatcher
# ---------------------------------------------------------------

def b2_probe(p, min_k, max_k):
    """Return (q, k) if some factor q = 2kp+1 with min_k ≤ k ≤ max_k
    divides M_p; else None."""
    if p in SMALL_KNOWN_P:
        return None
    return _py_b2(p, min_k, max_k)


def _py_b2(p, min_k, max_k):
    mp = (1 << p) - 1 if p < 64 else None
    two_p = 2 * p
    for k in range(max(min_k, 1), max_k + 1):
        q = two_p * k + 1
        if mp is not None and q >= mp:
            break
        if q % 24 not in (1, 7, 17, 23):
            continue
        if not is_prime_64(q):
            continue
        if pow(2, p, q) == 1:
            return (q, k)
    return None

# ---------------------------------------------------------------
# Initial filter
# ---------------------------------------------------------------

def initial_filter(start_p, end_p):
    out = []
    for p in range(start_p, end_p + 1):
        if p in KNOWN_COMPOSITE_DOUBLE_MERSENNE_EXPONENTS:
            continue
        if not gate0(p):
            continue
        if p in SMALL_KNOWN_P:
            out.append(p); continue
        if is_prime_64(p):
            out.append(p)
    return out

# ---------------------------------------------------------------
# Deepening loop
# ---------------------------------------------------------------

def deepen(survivors, start_k, multiplier, interactive, max_rounds=None):
    prev_k = 0
    k = start_k
    rnd = 0
    history = []

    while survivors and (max_rounds is None or rnd < max_rounds):
        rnd += 1
        t0 = time.time()

        new_survivors = []
        eliminated = []
        probed_actual = 0
        for p in survivors:
            if p in KNOWN_MERSENNE_EXPONENTS:
                new_survivors.append(p)
                continue
            probed_actual += 1
            hit = b2_probe(p, prev_k + 1, k)
            if hit is None:
                new_survivors.append(p)
            else:
                q, kk = hit
                eliminated.append((p, q, kk))

        elapsed = time.time() - t0
        killed = len(eliminated)
        print()
        print("-" * 72)
        print(f"Round {rnd}  |  k ∈ ({prev_k}, {k}]  |  "
              f"probed {probed_actual}  killed {killed}  "
              f"survive {len(new_survivors)}  |  {elapsed:.1f}s")
        if eliminated:
            ks = sorted(kk for _, _, kk in eliminated)
            print(f"  eliminated at k: "
                  f"min={ks[0]}  median={ks[len(ks)//2]}  max={ks[-1]}")
        print("-" * 72)

        history.append({
            "round": rnd, "k_lo": prev_k + 1, "k_hi": k,
            "probed": probed_actual, "killed": killed,
            "survived": len(new_survivors), "seconds": elapsed,
            "eliminated": eliminated,
        })

        if not eliminated and not interactive:
            print("No new eliminations. Stopping.")
            break

        survivors = new_survivors
        if not survivors:
            print("All non-known candidates eliminated.")
            break

        prev_k = k

        if interactive:
            print(f"\nSurvivors remaining: {len(survivors)}")
            try:
                ans = input(f"Next max_k [{k * multiplier}] ('q' to stop): ").strip()
            except EOFError:
                break
            if ans.lower() in ("q", "quit", "exit"):
                break
            k = int(ans) if ans else k * multiplier
        else:
            k *= multiplier

    return survivors, history

# ---------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------

def elimination_profile(history):
    all_elims = [el for h in history for el in h["eliminated"]]
    if not all_elims:
        return
    from collections import Counter
    buckets = Counter(int(math.log2(kk)) if kk > 0 else 0
                      for _, _, kk in all_elims)
    print("\nElimination profile (k bucket, log2):")
    for b in sorted(buckets):
        lo, hi = 1 << b, (1 << (b + 1)) - 1
        bar = "#" * min(60, buckets[b])
        print(f"  [{lo:>9}, {hi:>9}]  {buckets[b]:>5}  {bar}")


def report(survivors):
    print()
    print("=" * 72)
    print(f"FINAL SURVIVORS: {len(survivors)}")
    print("=" * 72)
    known = sorted(p for p in survivors if p in KNOWN_MERSENNE_EXPONENTS)
    rest  = sorted(p for p in survivors if p not in KNOWN_MERSENNE_EXPONENTS)
    if known:
        print(f"  Known Mersenne exponents captured ({len(known)}): {known}")
    if rest:
        print(f"  Unverified candidates ({len(rest)}):")
        for p in rest[:50]:
            print(f"    {p}")
        if len(rest) > 50:
            print(f"    ... and {len(rest) - 50} more")
    print("=" * 72)

# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Iterative deepening filter.")
    ap.add_argument("--start", type=int, default=101)
    ap.add_argument("--end",   type=int, default=500)
    ap.add_argument("--k",     type=int, default=1000,
                    help="initial max_k for round 1")
    ap.add_argument("--mult",  type=int, default=10,
                    help="k multiplier per round (non-interactive)")
    ap.add_argument("--rounds", type=int, default=None)
    ap.add_argument("--interactive", action="store_true",
                    help="prompt for max_k each round")
    ap.add_argument("--save", type=str, default=None)
    args = ap.parse_args()

    if False:  # Pure Python backend
        print(f"Backend: pure Python")
    else:
        print(f"B2 backend: python (slower)")

    print(f"Range: [{args.start}, {args.end}]  initial_k={args.k}  "
          f"multiplier={args.mult}")

    t0 = time.time()
    survivors = initial_filter(args.start, args.end)
    print(f"Initial filter: {len(survivors)} survivors in {time.time() - t0:.2f}s")

    survivors, history = deepen(
        survivors, start_k=args.k, multiplier=args.mult,
        interactive=args.interactive, max_rounds=args.rounds,
    )

    elimination_profile(history)
    report(survivors)

    if args.save:
        p = Path(args.save)
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "w") as f:
            json.dump({
                "start": args.start, "end": args.end,
                "initial_k": args.k, "history": history,
                "survivors": sorted(survivors),
                "backend": _RUST_SIG,
            }, f, indent=2)
        print(f"Saved -> {p}")


if __name__ == "__main__":
    main()