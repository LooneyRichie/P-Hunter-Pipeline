"""
p_hunter_eliminate.py — push farm_survivors.json deeper to eliminate more p, plus spacing stats.

    python p_hunter_eliminate.py run                    # 60 CPU-min cap per p, all cores
    python p_hunter_eliminate.py run --max-minutes 30 --max-hours 8
    python p_hunter_eliminate.py spacing                # where the next Mersenne exponent likely is

`run` reads farm_survivors.json (written by p_hunter_farm.py), trial-factors every survivor in
common-depth rounds (so budget is spread evenly), and stops a p once it has used its per-p
CPU budget. Progress lives in eliminate_state.json; results go to farm_survivors_deepened.*.
Survivors stay unproven: only Lucas-Lehmer / PRP can prove M_p prime.
"""
import argparse
import concurrent.futures as cf
import json
import math
import multiprocessing as mp
import os
import sys
import time

import revised_binV3 as core

try:
    from p_hunter_tf import tf_window
except ImportError:
    tf_window = None

SRC_JSON = "farm_survivors.json"
STATE_FILE = "eliminate_state.json"
OUT_JSON = "farm_survivors_deepened.json"
OUT_TXT = "farm_survivors_deepened.txt"
BAND = 10**9                 # k per task, about 40 s at the measured rate
MIN_STEP = 5 * 10**7         # smaller remaining gaps are not worth a round
DEFAULT_RATE = 4e-8          # seconds per k per core, from the p=136279901 benchmark
GROWTH = 3                   # common depth multiplier per round
EGAMMA = 1.7810724179901979  # e^gamma
LAMBDA = EGAMMA / math.log(2)  # Wagstaff: expected Mersenne exponents <= x is LAMBDA * ln x


# ---------------------------------------------------------------- state

def save_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1)
    os.replace(tmp, path)


def merge(iv):
    out = []
    for lo, hi in sorted(iv):
        if out and lo <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return out


def depth(entry):
    cov = merge(entry["covered"])
    return cov[0][1] if cov and cov[0][0] == 1 else 0


def load_state():
    if os.path.exists(STATE_FILE):
        with open(STATE_FILE) as f:
            return json.load(f)
    return {"rate": DEFAULT_RATE, "entries": {}}


def sync_from_farm(state):
    if not os.path.exists(SRC_JSON):
        sys.exit(f"{SRC_JSON} not found; run p_hunter_farm.py first")
    with open(SRC_JSON) as f:
        src = json.load(f)
    added = 0
    for s in src["survivors"]:
        key = str(s["p"])
        if key not in state["entries"]:
            state["entries"][key] = {"covered": [[1, s["depth_k"]]], "spent": 0.0,
                                     "status": "alive", "factor": None}
            added += 1
    return added


def write_outputs(state):
    alive, dead = [], []
    for key, e in state["entries"].items():
        p = int(key)
        d = depth(e)
        row = {"p": p, "depth_k": d, "bits": round(math.log2(2 * max(d, 1) * p + 1), 2),
               "score": core.mersenne_likeness_score(p)[0],
               "cpu_min": round(e["spent"] / 60, 1), "status": e["status"]}
        if e["status"] == "dead":
            row["factor"] = e["factor"]
            dead.append(row)
        else:
            alive.append(row)
    alive.sort(key=lambda r: (-r["bits"], -r["score"], r["p"]))
    dead.sort(key=lambda r: r["p"])
    save_json(OUT_JSON, {"alive": alive, "eliminated": dead})
    with open(OUT_TXT, "w") as f:
        f.write("# alive = no factor up to depth_k (unproven). status 'capped' = per-p budget used.\n")
        f.write(f"# {'p':>12} {'depth_k':>14} {'bits':>6} {'score':>5} {'cpu_min':>8} status\n")
        for r in alive:
            f.write(f"  {r['p']:>12} {r['depth_k']:>14} {r['bits']:>6.1f} {r['score']:>5} "
                    f"{r['cpu_min']:>8.1f} {r['status']}\n")
        f.write(f"\n# eliminated ({len(dead)})\n")
        for r in dead:
            f.write(f"  {r['p']:>12}  factor {r['factor'][0]}  (k={r['factor'][1]})\n")
    return alive, dead


# ---------------------------------------------------------------- run

def _band(args):
    p, lo, hi = args
    t0 = time.perf_counter()
    hit = tf_window(p, lo, hi)
    return p, lo, hi, hit, time.perf_counter() - t0


def plan_round(state, cap_s):
    """Return (target_depth, tasks, capped_keys) for the next common-depth round, or None."""
    ents = state["entries"]
    rate = state["rate"]
    live = [(k, e) for k, e in ents.items() if e["status"] == "alive"]
    if not live:
        return None
    floor = min(depth(e) for _, e in live)
    target = max(floor * GROWTH, floor + MIN_STEP)
    tasks, capped = [], []
    for key, e in live:
        p, d = int(key), depth(e)
        budget_k = int(max(cap_s - e["spent"], 0) / rate)
        t = min(target, d + budget_k)
        if t - d < MIN_STEP:
            if d < target:
                capped.append(key)
            continue
        lo = d + 1
        while lo <= t:
            hi = min(lo + BAND - 1, t)
            tasks.append((p, lo, hi))
            lo = hi + 1
    tasks.sort(key=lambda x: (x[1], x[0]))
    return target, tasks, capped


def run(args):
    if tf_window is None:
        sys.exit("p_hunter_tf missing: cd p_hunter_tf && maturin develop --release")
    state = load_state()
    added = sync_from_farm(state)
    save_json(STATE_FILE, state)
    cap_s = args.max_minutes * 60
    deadline = time.time() + args.max_hours * 3600 if args.max_hours else float("inf")
    ents = state["entries"]
    n_alive = sum(e["status"] == "alive" for e in ents.values())
    print(f"workers={args.workers}  cap={args.max_minutes} CPU-min/p  "
          f"alive={n_alive} (+{added} new)  rate={state['rate']:.2e} s/k")

    try:
        ctx = mp.get_context("fork")
    except ValueError:
        ctx = mp.get_context()
    pool = cf.ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx)
    rnd = 0
    try:
        while time.time() < deadline:
            plan = plan_round(state, cap_s)
            if plan is None:
                print("nothing alive.")
                break
            target, tasks, capped = plan
            for key in capped:
                ents[key]["status"] = "capped"
            if capped:
                print(f"  {len(capped)} p reached the {args.max_minutes:g}-min cap")
            if not tasks:
                break
            rnd += 1
            alive_before = sum(e["status"] == "alive" for e in ents.values())
            est = len(tasks) * BAND * state["rate"] / args.workers
            print(f"\nRound {rnd}: depth -> k<={target:,}  p={len({t[0] for t in tasks})}  "
                  f"est wall {est / 60:.0f} min")
            t0, cpu, killed = time.time(), 0.0, 0
            it, pending, stop = iter(tasks), {}, False
            while True:
                while len(pending) < args.workers * 2 and not stop:
                    t = next(it, None)
                    if t is None:
                        break
                    if ents[str(t[0])]["status"] != "alive":
                        continue
                    if time.time() > deadline:
                        stop = True
                        break
                    pending[pool.submit(_band, t)] = t
                if not pending:
                    break
                done, _ = cf.wait(pending, return_when=cf.FIRST_COMPLETED)
                for f in done:
                    pending.pop(f)
                    p, lo, hi, hit, dt = f.result()
                    e = ents[str(p)]
                    if e["status"] != "alive":
                        continue
                    e["spent"] += dt
                    cpu += dt
                    if hit is not None:
                        q, k = hit
                        if q == 2 * k * p + 1 and pow(2, p, q) == 1:
                            e["status"], e["factor"] = "dead", [q, k]
                            killed += 1
                            print(f"  [-] p={p}  factor {q}  (k={k:,})")
                        else:
                            print(f"  [!] p={p} backend returned bad factor {hit}; ignored")
                    else:
                        e["covered"] = merge(e["covered"] + [[lo, hi]])
                        state["rate"] = 0.8 * state["rate"] + 0.2 * dt / (hi - lo + 1)
                save_json(STATE_FILE, state)
            write_outputs(state)
            wall = time.time() - t0
            print(f"  round {rnd}: killed {killed}/{alive_before}  wall {wall / 60:.1f} min  "
                  f"cpu {cpu / 3600:.2f} h  kills/cpu-hour {killed / max(cpu / 3600, 1e-9):.2f}")
            if stop:
                print("deadline reached.")
                break
    except KeyboardInterrupt:
        print("\ninterrupted; completed bands are saved.")
        pool.shutdown(wait=False, cancel_futures=True)
    else:
        pool.shutdown()
    alive, dead = write_outputs(state)
    save_json(STATE_FILE, state)
    print(f"\nalive={len(alive)}  eliminated={len(dead)}  -> {OUT_TXT}")


# ---------------------------------------------------------------- spacing

def spacing(args):
    ks = sorted(core.KNOWN_MERSENNE_EXPONENTS)
    p0 = args.after or ks[-1]
    print(f"Known Mersenne exponents: {len(ks)}  (largest {ks[-1]:,})")
    print(f"Wagstaff density: {LAMBDA:.3f} exponents per unit of ln p "
          f"(predicts {LAMBDA * math.log(ks[-1]):.1f} up to {ks[-1]:,}; actual {len(ks)})\n")

    print("Recent gaps (ratio = p_n / p_n-1):")
    for a, b in list(zip(ks, ks[1:]))[-12:]:
        print(f"  {a:>11,} -> {b:>11,}   ratio {b / a:5.3f}   gap {b - a:>11,}")
    lr = sorted(math.log(b / a) for a, b in zip(ks, ks[1:]) if a >= 1000)
    print(f"\nln-ratio since p>=1000: mean {sum(lr) / len(lr):.3f} (theory 1/lambda = {1 / LAMBDA:.3f}), "
          f"median {lr[len(lr) // 2]:.3f} (theory {math.log(2) / LAMBDA:.3f})")

    print(f"\nNext exponent after p0={p0:,} (memoryless in ln p; skipping a few thousand changes nothing):")
    print(f"  {'cum. prob':>9}  p <=")
    for q in (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95):
        r = math.exp(-math.log(1 - q) / LAMBDA)
        print(f"  {q:>9.0%}  {p0 * r:>14,.0f}   (x{r:.3f})")

    print("\nBands: expected exponents in band, primes to test, chance per prime tested")
    print(f"  {'band':>25}  {'E[count]':>9}  {'primes':>11}  {'P(prime) each':>13}")
    edges = [p0, 1.5e8, 2e8, 3e8, 5e8, 1e9]
    for a, b in zip(edges, edges[1:]):
        if b <= a:
            continue
        exp_count = LAMBDA * math.log(b / a)
        n_primes = (b - a) / math.log(math.sqrt(a * b))
        print(f"  {a:>11,.0f}..{b:>11,.0f}  {exp_count:>9.3f}  {n_primes:>11,.0f}  "
              f"{exp_count / n_primes:>13.2e}")
    print("  Per-prime odds fall as 1/p, so the best return per test is the lowest untested p;")
    print("  no band is more concentrated than another.")

    print("\nTrial factoring to q <= B raises the odds by about ln(B)/ln(2p) per survivor:")
    for kmax, lab in ((10**9, "k<=1e9"), (9 * 10**10, "k<=9e10 (~1 CPU-hour/p)"),
                      (10**13, "k<=1e13"), (10**16, "k<=1e16 (GPU territory)")):
        qb = 2 * kmax * p0
        surv = math.log(2 * p0) / math.log(qb)
        print(f"  {lab:<26} q~2^{math.log2(qb):5.1f}  survivors {surv:6.1%}  odds boost x{1 / surv:4.2f}")
    per_prime = LAMBDA * math.log(p0) / p0
    print(f"\nBase odds that one prime p near {p0:,} is a Mersenne exponent: {per_prime:.2e}")
    print(f"Testing 100 such survivors (boost ~x2.5): about {100 * per_prime * 2.5:.1e} total chance.")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="deepen trial factoring on farm survivors")
    r.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    r.add_argument("--max-minutes", type=float, default=60, help="CPU-minute cap per p")
    r.add_argument("--max-hours", type=float, default=None, help="wall-clock limit for this run")
    s = sub.add_parser("spacing", help="Mersenne exponent spacing / next-exponent estimate")
    s.add_argument("--after", type=int, default=None, help="p0 (default: largest known exponent)")
    args = ap.parse_args()
    (run if args.cmd == "run" else spacing)(args)


if __name__ == "__main__":
    main()
