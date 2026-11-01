"""
p_hunter_farm.py — multi-worker survivor farm for Mersenne exponents above 136279841.

Walks prime p upward, runs staged trial factoring (q = 2kp+1) on each p across
worker processes, and stops once TARGET new survivors (no factor up to depth_k)
exist. State is checkpointed so the next run continues for TARGET more.

    python p_hunter_farm.py                       # first 100, depth k<=1e9
    python p_hunter_farm.py --resume              # next 100 from checkpoint
    python p_hunter_farm.py --depth-k 10000000000 --workers 4

Survivors are NOT proven prime: only Lucas-Lehmer / PRP (GIMPS software) can do that.
Likeness score is ranking only, never a filter.
"""
import argparse
import json
import math
import multiprocessing as mp
import os
import sys
import time
import concurrent.futures as cf

import revised_binV3 as core

try:
    from p_hunter_tf import tf_window as _tf_window
    BACKEND = "p_hunter_tf"
except ImportError:
    _tf_window = None
    try:
        from p_hunter_rs import gate_b2_factor_probe as _rs_window
        BACKEND = "p_hunter_rs"
    except ImportError:
        _rs_window = None
        BACKEND = "none"

DEFAULT_START = 136279843
STATE_FILE = "farm_state.json"
OUT_JSON = "farm_survivors.json"
OUT_TXT = "farm_survivors.txt"
LADDER = (10**6, 10**8, 10**9, 10**10, 6 * 10**10)
RS_MAX_K = 60_000_000_000  # p_hunter_rs skips q > 2^64 silently


def window(p, k_lo, k_hi):
    """Return (q, k) if some q = 2kp+1 with k in [k_lo, k_hi] divides M_p, else None."""
    if _tf_window is not None:
        return _tf_window(p, k_lo, k_hi)
    return _rs_window(p, k_lo, k_hi)


def stages_for(depth_k):
    return [k for k in LADDER if k < depth_k] + [depth_k]


def probe_p(args):
    """Worker: climb the ladder for one p. Returns (p, factor_or_None)."""
    p, depth_k = args
    lo = 1
    for hi in stages_for(depth_k):
        hit = window(p, lo, hi)
        if hit is not None:
            return p, hit
        lo = hi + 1
    return p, None


def depth_bits(p, depth_k):
    return math.log2(2 * depth_k * p + 1)


def next_primes(start, count):
    """Next `count` primes >= start that pass the wheel; returns (primes, next_unscanned)."""
    out = []
    p = start
    while len(out) < count:
        if core.gate0_candidate_stream(p) and core.is_prime(p):
            out.append(p)
        p += 1
    return out, p


def load_state(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def save_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1)
    os.replace(tmp, path)


def rank(survivors):
    """Deepest first, then likeness score, then p."""
    for s in survivors:
        s["score"], _ = core.mersenne_likeness_score(s["p"])
    return sorted(survivors, key=lambda s: (-s["bits"], -s["score"], s["p"]))


def write_outputs(state):
    ranked = rank(state["survivors"])
    start_p = state.get("start_p")
    swept = [[start_p, state["next_p"] - 1]] if start_p is not None else []
    save_json(OUT_JSON, {"next_p": state["next_p"],
                         "prior_ranges_swept": swept, "survivors": ranked})
    with open(OUT_TXT, "w") as f:
        f.write("# unproven survivors: no factor q=2kp+1 with k <= depth_k. Test with LL/PRP.\n")
        f.write(f"# {'p':>12}  {'depth_k':>13}  {'bits':>5}  score  binary\n")
        for s in ranked:
            f.write(f"  {s['p']:>12}  {s['depth_k']:>13}  {s['bits']:5.1f}  {s['score']:>5}  {bin(s['p'])[2:]}\n")
    return ranked


def make_pool(workers):
    try:
        ctx = mp.get_context("fork")
    except ValueError:
        ctx = mp.get_context()
    return cf.ProcessPoolExecutor(max_workers=workers, mp_context=ctx)


def main():
    ap = argparse.ArgumentParser(description="Farm trial-factoring survivors above 136279841.")
    ap.add_argument("--start", type=int, default=DEFAULT_START)
    ap.add_argument("--target", type=int, default=100, help="new survivors to find this run")
    ap.add_argument("--depth-k", type=int, default=10**9, help="trial-factor every k <= this")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--chunk", type=int, default=None, help="primes per batch (default 4*workers)")
    ap.add_argument("--resume", action="store_true", help="continue from the checkpoint")
    ap.add_argument("--state", default=STATE_FILE)
    args = ap.parse_args()

    if BACKEND == "none":
        sys.exit("No Rust backend. Build p_hunter_tf: cd p_hunter_tf && maturin develop --release")
    if BACKEND == "p_hunter_rs":
        print("[warn] p_hunter_tf missing; using p_hunter_rs (q capped at 2^64, k_max bytes per worker).")
        if args.depth_k > RS_MAX_K:
            sys.exit(f"--depth-k above {RS_MAX_K:,} is not rigorous with p_hunter_rs")
    if args.start <= 136279841:
        sys.exit("--start must be above 136279841")

    state = load_state(args.state) if args.resume else None
    if args.resume and state is None:
        sys.exit(f"--resume given but no checkpoint at {args.state}")
    if state is None:
        if os.path.exists(args.state):
            sys.exit(f"{args.state} exists; use --resume or remove it to start over")
        state = {"start_p": args.start, "next_p": args.start, "survivors": [],
                 "eliminated": 0, "primes_scanned": 0}
    else:
        state.setdefault("start_p", None)

    chunk = args.chunk or 4 * args.workers
    print(f"backend={BACKEND}  workers={args.workers}  depth_k={args.depth_k:,} "
          f"(~2^{depth_bits(state['next_p'], args.depth_k):.1f})  target={args.target}")
    print(f"start p={state['next_p']:,}  stored survivors={len(state['survivors'])}")

    found = []
    t0 = time.time()
    pool = make_pool(args.workers)
    try:
        while len(found) < args.target:
            primes, scanned_to = next_primes(state["next_p"], chunk)
            tc = time.time()
            results = list(pool.map(probe_p, [(p, args.depth_k) for p in primes]))
            alive = sorted(p for p, hit in results if hit is None)
            killed = len(primes) - len(alive)

            need = args.target - len(found)
            if len(alive) >= need:
                # Cut at the need-th survivor so resume continues exactly after it.
                alive = alive[:need]
                cut = alive[-1]
                kept = [p for p in primes if p <= cut]
                killed = len(kept) - len(alive)
                state["next_p"] = cut + 1
                primes = kept
            else:
                state["next_p"] = scanned_to

            for p in alive:
                rec = {"p": p, "depth_k": args.depth_k, "bits": round(depth_bits(p, args.depth_k), 2)}
                found.append(rec)
                state["survivors"].append(rec)
            state["eliminated"] += killed
            state["primes_scanned"] += len(primes)
            save_json(args.state, state)
            write_outputs(state)
            print(f"  p<{state['next_p']:,}  tested={len(primes)} killed={killed} "
                  f"survived={len(alive)}  run total={len(found)}/{args.target}  "
                  f"({time.time() - tc:.1f}s, {time.time() - t0:.0f}s elapsed)")
    except KeyboardInterrupt:
        print("\ninterrupted; checkpoint holds the last finished batch.")
        pool.shutdown(wait=False, cancel_futures=True)
        return
    pool.shutdown()

    ranked = write_outputs(state)
    print()
    print("=" * 72)
    print(f"{len(found)} new survivors; {len(state['survivors'])} total; next_p={state['next_p']:,}")
    print(f"scanned {state['primes_scanned']} primes, eliminated {state['eliminated']}")
    print("Top 20 (depth, then likeness):")
    for s in ranked[:20]:
        print(f"  p={s['p']:<12} 2^{s['bits']:.1f}  score={s['score']}/10")
    print(f"Full list: {OUT_TXT} / {OUT_JSON}.  Resume for 100 more: python p_hunter_farm.py --resume")
    print("Survivors are unproven; run LL/PRP (e.g. GIMPS software) on the top entries.")
    print("=" * 72)


if __name__ == "__main__":
    main()
