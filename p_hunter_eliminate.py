"""
p_hunter_eliminate.py — menu-driven, one-p-at-a-time deep trial factoring of farm survivors.

Run:  python p_hunter_eliminate.py

Works down the list from farm_survivors.json. Each p is trial-factored (q = 2kp+1) with all
workers on that single p, at max k = START, then doubled (START, 2x, 4x, ...) until it is
eliminated, hits the k ceiling, or uses up its wall-clock time limit. Eliminated p are crossed
off; the rest are marked "passed" for further testing. State is saved after every band and
after every p, so Ctrl+C is always safe.

Survivors stay unproven: only Lucas-Lehmer / PRP can prove M_p prime.
"""
import concurrent.futures as cf
import json
import math
import multiprocessing as mp
import os
import signal
import sys
import time
from datetime import datetime, timezone

import revised_binV3 as core

try:
    from p_hunter_tf import tf_window
except ImportError:
    tf_window = None

SRC_JSON = "farm_survivors.json"
STATE_FILE = "eliminate_state.json"
OUT_JSON = "farm_survivors_deepened.json"
OUT_TXT = "farm_survivors_deepened.txt"

BAND = 10**9                 # k per task, about 40-60 s of one core
MIN_BAND = 2 * 10**7         # smaller bands are not worth dispatching
DEFAULT_RATE = 5e-8          # seconds per k per core until measured
MAX_CEILING = 10**18         # keeps k inside u64 and q below 2^100 for p < 1e9
HEARTBEAT = 120              # seconds between progress lines inside a level
EGAMMA = 1.7810724179901979
LAMBDA = EGAMMA / math.log(2)  # Wagstaff: expected Mersenne exponents <= x is LAMBDA * ln x

DEFAULTS = {
    "minutes": 60.0,
    "start_k": 50_000_000_000,
    "ceiling_k": 10**14,
    "workers": os.cpu_count() or 1,
}


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


def bits(p, k):
    return math.log2(2 * max(k, 1) * p + 1)


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def load_state():
    state = {"rate": DEFAULT_RATE, "entries": {}}
    if os.path.exists(STATE_FILE):
        with open(STATE_FILE) as f:
            state = json.load(f)
    state["settings"] = {**DEFAULTS, **state.get("settings", {})}
    legacy = {"alive": "pending", "capped": "passed"}   # older eliminator wording
    for e in state["entries"].values():
        e["status"] = legacy.get(e["status"], e["status"])
    return state


def sync_from_farm(state):
    """Append survivors that appeared in farm_survivors.json since last time."""
    if not os.path.exists(SRC_JSON):
        return 0
    with open(SRC_JSON) as f:
        src = json.load(f)
    swept = src.get("prior_ranges_swept")
    if swept is None and os.path.exists("farm_state.json"):
        with open("farm_state.json") as f:
            farm_state = json.load(f)
        start_p = farm_state.get("start_p")
        next_p = farm_state.get("next_p")
        swept = [[start_p, next_p - 1]] if start_p is not None and next_p is not None else []
    state["prior_ranges_swept"] = swept or state.get("prior_ranges_swept", [])
    added = 0
    for s in src["survivors"]:
        key = str(s["p"])
        if key not in state["entries"]:
            state["entries"][key] = {"covered": [[1, s["depth_k"]]], "spent": 0.0,
                                     "status": "pending", "factor": None}
            added += 1
    return added


def counts(state):
    c = {"pending": 0, "dead": 0, "passed": 0}
    for e in state["entries"].values():
        c[e["status"]] += 1
    return c


def write_outputs(state):
    rows = {"pending": [], "passed": [], "dead": []}
    tool = "p_hunter_tf 0.1.0"
    method = "2kp+1 sieve, q=2kp+1"
    for key, e in state["entries"].items():
        p = int(key)
        d = depth(e)
        sieve_bits = (2 * d * p + 1).bit_length()
        if e["status"] == "dead":
            notes = f"Eliminated by factor q={e['factor'][0]} at k={e['factor'][1]}."
        elif e["status"] == "passed":
            notes = f"No factor found through k={d:,} (below 2^{sieve_bits}); M_p remains unproven."
        else:
            notes = f"Current contiguous coverage through k={d:,}; candidate not yet fully processed."
        rows[e["status"]].append({
            "p": p, "sieve_bits": sieve_bits, "k_max": d,
            "tool": tool, "method": method,
            "timestamp": e.get("updated_at"),
            "prior_ranges_swept": state.get("prior_ranges_swept", []),
            "notes": notes,
            "depth_k": d, "bits": round(bits(p, d), 2),
            "score": core.mersenne_likeness_score(p)[0],
            "minutes": round(e["spent"] / 60, 1), "factor": e["factor"]})
    for lst in (rows["passed"], rows["pending"]):
        lst.sort(key=lambda r: (-r["bits"], -r["score"], r["p"]))
    rows["dead"].sort(key=lambda r: r["p"])
    save_json(OUT_JSON, rows)
    with open(OUT_TXT, "w") as f:
        f.write("# PASSED = no factor up to depth_k (unproven). Test these with LL/PRP or GPU factoring.\n")
        f.write(f"# {'p':>12} {'sieve_bits':>10} {'k_max':>16} {'score':>5} {'minutes':>8}\n")
        for tag in ("passed", "pending"):
            f.write(f"\n## {tag.upper()} ({len(rows[tag])})\n")
            for r in rows[tag]:
                f.write(f"  {r['p']:>12} {r['sieve_bits']:>10} {r['k_max']:>16} {r['score']:>5} {r['minutes']:>8.1f}\n")
        f.write(f"\n## ELIMINATED ({len(rows['dead'])})\n")
        for r in rows["dead"]:
            f.write(f"  {r['p']:>12}  factor {r['factor'][0]}  (k={r['factor'][1]:,})\n")


# ---------------------------------------------------------------- factoring

def _worker_init():
    signal.signal(signal.SIGINT, signal.SIG_IGN)   # main process handles Ctrl+C


def _band(args):
    p, lo, hi = args
    t0 = time.perf_counter()
    hit = tf_window(p, lo, hi)
    return p, lo, hi, hit, time.perf_counter() - t0


def make_pool(workers):
    try:
        ctx = mp.get_context("fork")
    except ValueError:
        ctx = mp.get_context()
    return cf.ProcessPoolExecutor(max_workers=workers, mp_context=ctx, initializer=_worker_init)


def levels(start, ceiling):
    out, k = [], start
    while k < ceiling:
        out.append(k)
        k *= 2
    out.append(ceiling)
    return out


def fmt_time(s):
    s = int(s)
    return f"{s // 3600}h{s % 3600 // 60:02d}m" if s >= 3600 else f"{s // 60}m{s % 60:02d}s"


def process_p(pool, state, key):
    """Deepen one p until eliminated, ceiling reached, or time limit used. Returns final status."""
    cfg = state["settings"]
    e = state["entries"][key]
    p = int(key)
    limit_s = cfg["minutes"] * 60
    base_spent = e["spent"]
    t_start = time.time()
    workers = cfg["workers"]
    last_beat = t_start
    out_of_time = False

    def elapsed_total():
        return base_spent + (time.time() - t_start)

    for lvl in levels(cfg["start_k"], cfg["ceiling_k"]):
        if depth(e) >= lvl:
            continue
        if elapsed_total() >= limit_s:
            break
        print(f"    k <= {lvl:>18,}  (~2^{bits(p, lvl):.1f})   spent {fmt_time(elapsed_total())}"
              f" / {fmt_time(limit_s)}")
        next_lo = depth(e) + 1
        pending = {}
        while True:
            while len(pending) < workers and next_lo <= lvl and not out_of_time:
                remaining = limit_s - elapsed_total()
                width = min(BAND, lvl - next_lo + 1, int(remaining / state["rate"]))
                if width < min(MIN_BAND, lvl - next_lo + 1):
                    out_of_time = True
                    break
                pending[pool.submit(_band, (p, next_lo, next_lo + width - 1))] = None
                next_lo += width
            if not pending:
                break
            done, _ = cf.wait(pending, return_when=cf.FIRST_COMPLETED)
            for f in done:
                del pending[f]
                _, lo, hi, hit, dt = f.result()
                if hit is not None:
                    q, k = hit
                    if q == 2 * k * p + 1 and pow(2, p, q) == 1:
                        for other in pending:
                            other.cancel()
                        e["status"], e["factor"] = "dead", [q, k]
                        e["spent"] = elapsed_total()
                        e["updated_at"] = utc_now()
                        return "dead"
                    print(f"    [!] backend returned a bad factor {hit}; ignoring it")
                    continue
                e["covered"] = merge(e["covered"] + [[lo, hi]])
                e["updated_at"] = utc_now()
                state["rate"] = 0.8 * state["rate"] + 0.2 * dt / (hi - lo + 1)
            e["spent"] = elapsed_total()
            save_json(STATE_FILE, state)
            if time.time() - last_beat >= HEARTBEAT:
                last_beat = time.time()
                print(f"      ... clear to k={depth(e):,}  (~2^{bits(p, depth(e)):.1f})  "
                      f"spent {fmt_time(elapsed_total())}")
        if out_of_time:
            break

    e["spent"] = elapsed_total()
    e["status"] = "passed"
    e["updated_at"] = utc_now()
    return "passed"


def run_eliminator(state):
    if tf_window is None:
        print("p_hunter_tf is not installed: cd p_hunter_tf && maturin develop --release")
        return
    added = sync_from_farm(state)
    if added:
        print(f"Picked up {added} new survivor(s) from {SRC_JSON}.")
    todo = [k for k, e in state["entries"].items() if e["status"] == "pending"]
    if not todo:
        print("Nothing pending. Use 'Re-queue' or run the farm for more survivors.")
        return
    cfg = state["settings"]
    print(f"\n{len(todo)} pending. Limit {cfg['minutes']:g} min/p, doubling from k={cfg['start_k']:,} "
          f"to ceiling {cfg['ceiling_k']:,}, {cfg['workers']} workers. Ctrl+C stops safely.\n")
    pool = make_pool(cfg["workers"])
    try:
        for i, key in enumerate(todo, 1):
            p = int(key)
            e = state["entries"][key]
            print(f"[{i}/{len(todo)}] p={p:,}  already clear to k={depth(e):,}")
            status = process_p(pool, state, key)
            save_json(STATE_FILE, state)
            write_outputs(state)
            if status == "dead":
                q, k = e["factor"]
                print(f"  X  ELIMINATED  factor q={q:,}  (k={k:,}, ~2^{math.log2(q):.1f})")
            else:
                print(f"  OK PASSED      clear to k={depth(e):,}  (~2^{bits(p, depth(e)):.1f})"
                      f"  in {fmt_time(e['spent'])}")
        print("\nFinished the list.")
    except KeyboardInterrupt:
        save_json(STATE_FILE, state)
        write_outputs(state)
        print("\nStopped. Finished bands are saved; the current p resumes where it left off.")
        pool.shutdown(wait=False, cancel_futures=True)
        return
    pool.shutdown()
    c = counts(state)
    print(f"Eliminated {c['dead']}  Passed {c['passed']}  Pending {c['pending']}  -> {OUT_TXT}")


# ---------------------------------------------------------------- menu helpers

def num(text):
    return int(float(text.replace(",", "").replace("_", "")))


def ask(prompt, default):
    raw = input(f"{prompt} [{default:,}]: " if isinstance(default, int) and default > 9999
                else f"{prompt} [{default}]: ").strip()
    return default if raw == "" else raw


def show_list(state):
    mark = {"pending": "[ ]", "passed": "[OK]", "dead": "[X]"}
    print()
    for key, e in state["entries"].items():
        p = int(key)
        d = depth(e)
        extra = ""
        if e["status"] == "dead":
            extra = f"factor {e['factor'][0]}"
        elif d:
            extra = f"clear to k={d:,} (~2^{bits(p, d):.1f})  {fmt_time(e['spent'])}"
        print(f"  {mark[e['status']]:<4} {p:>12}  {extra}")
    c = counts(state)
    print(f"\n  [X] eliminated {c['dead']}   [OK] passed {c['passed']}   [ ] pending {c['pending']}")


def edit_settings(state):
    cfg = state["settings"]
    print("\nSettings (Enter keeps the current value; numbers like 5e10 or 50,000,000,000 work)")
    try:
        minutes = float(ask("Time limit per p, minutes", cfg["minutes"]))
        start_k = num(str(ask("Starting max k", cfg["start_k"])))
        ceiling = num(str(ask("Ceiling max k", cfg["ceiling_k"])))
        workers = int(ask("Workers", cfg["workers"]))
    except ValueError:
        print("Not a number; settings unchanged.")
        return
    if minutes <= 0 or start_k < 1 or workers < 1 or not start_k <= ceiling <= MAX_CEILING:
        print(f"Invalid (need start <= ceiling <= {MAX_CEILING:,}); settings unchanged.")
        return
    cfg.update(minutes=minutes, start_k=start_k, ceiling_k=ceiling, workers=workers)
    save_json(STATE_FILE, state)
    lvls = levels(start_k, ceiling)
    print(f"Saved. Levels: {', '.join(f'{k:.3g}' for k in lvls[:8])}{' ...' if len(lvls) > 8 else ''}")
    rate = state["rate"]
    reach = rate and minutes * 60 * workers / rate
    print(f"At the measured {rate:.2e} s/k/core, {minutes:g} min x {workers} workers "
          f"reaches about k={reach:,.0f} per p.")


def requeue(state):
    passed = [k for k, e in state["entries"].items() if e["status"] == "passed"]
    if not passed:
        print("No passed p to re-queue.")
        return
    print(f"{len(passed)} passed p. Time is counted cumulatively per p, so raise the time limit")
    print("(Settings) first if you want them to get more time.")
    if input("Re-queue all passed p? [y/N]: ").strip().lower() == "y":
        for k in passed:
            state["entries"][k]["status"] = "pending"
        save_json(STATE_FILE, state)
        print(f"{len(passed)} re-queued.")


def spacing():
    ks = sorted(core.KNOWN_MERSENNE_EXPONENTS)
    raw = input(f"Estimate the next exponent after p0 [{ks[-1]:,}]: ").strip()
    try:
        p0 = num(raw) if raw else ks[-1]
    except ValueError:
        p0 = ks[-1]
    print(f"\nKnown Mersenne exponents: {len(ks)}  (largest {ks[-1]:,})")
    print(f"Wagstaff density: {LAMBDA:.3f} exponents per unit of ln p "
          f"(predicts {LAMBDA * math.log(ks[-1]):.1f} up to {ks[-1]:,}; actual {len(ks)})\n")
    print("Recent gaps:")
    for a, b in list(zip(ks, ks[1:]))[-10:]:
        print(f"  {a:>11,} -> {b:>11,}   ratio {b / a:5.3f}   gap {b - a:>11,}")
    lr = sorted(math.log(b / a) for a, b in zip(ks, ks[1:]) if a >= 1000)
    print(f"\nln-ratio since p>=1000: mean {sum(lr) / len(lr):.3f} (theory {1 / LAMBDA:.3f}), "
          f"median {lr[len(lr) // 2]:.3f} (theory {math.log(2) / LAMBDA:.3f})")
    print(f"\nNext exponent after {p0:,} (memoryless in ln p, so skipping a few thousand changes nothing):")
    for q in (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95):
        r = math.exp(-math.log(1 - q) / LAMBDA)
        print(f"  {q:>4.0%} chance it is at or below {p0 * r:>14,.0f}   (x{r:.3f})")
    print("\nBand: expected exponents, primes to test, chance per prime tested")
    edges = [p0, 1.5e8, 2e8, 3e8, 5e8, 1e9]
    for a, b in zip(edges, edges[1:]):
        if b <= a:
            continue
        ex = LAMBDA * math.log(b / a)
        n = (b - a) / math.log(math.sqrt(a * b))
        print(f"  {a:>13,.0f}..{b:>13,.0f}  {ex:>6.3f}  {n:>11,.0f}  {ex / n:.2e}")
    print("  Per-prime odds fall as 1/p: no band is richer, the lowest untested p is the best bet.")
    print("\nOdds boost from trial factoring to q <= B is about ln(B)/ln(2p) per survivor:")
    for kmax in (10**9, 5 * 10**10, 4 * 10**11, 10**13, 10**16):
        qb = 2 * kmax * p0
        s = math.log(2 * p0) / math.log(qb)
        print(f"  k<={kmax:<18,} q~2^{math.log2(qb):5.1f}  survive {s:6.1%}  boost x{1 / s:4.2f}")
    pp = LAMBDA * math.log(p0) / p0
    print(f"\nBase odds one prime near {p0:,} is a Mersenne exponent: {pp:.2e}")


# ---------------------------------------------------------------- menu

def main():
    state = load_state()
    if sync_from_farm(state):
        save_json(STATE_FILE, state)
    while True:
        c = counts(state)
        cfg = state["settings"]
        print()
        print("=" * 70)
        print("P_HUNTER ELIMINATOR")
        print("=" * 70)
        print(f"List     : {len(state['entries'])} survivors   "
              f"eliminated {c['dead']}   passed {c['passed']}   pending {c['pending']}")
        print(f"Settings : {cfg['minutes']:g} min/p, k {cfg['start_k']:,} doubling to "
              f"{cfg['ceiling_k']:,}, {cfg['workers']} workers")
        print("-" * 70)
        print("1. Start / continue eliminating (works down the list)")
        print("2. Show the list")
        print("3. Settings")
        print("4. Re-queue passed p for deeper testing")
        print("5. Where might the next exponent be? (spacing)")
        print("0. Exit")
        print("=" * 70)
        try:
            choice = input("Choice [0-5]: ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            break
        try:
            if choice == "1":
                run_eliminator(state)
            elif choice == "2":
                show_list(state)
            elif choice == "3":
                edit_settings(state)
            elif choice == "4":
                requeue(state)
            elif choice == "5":
                spacing()
            elif choice == "0":
                break
            else:
                print("invalid.")
        except KeyboardInterrupt:
            print()


if __name__ == "__main__":
    main()
