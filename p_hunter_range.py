"""
p_hunter_range.py — scan a chosen p range with p_hunter_farm's factoring and collect survivors.

Run:  python p_hunter_range.py              (prompts for start p, end p, depth, ...)
      python p_hunter_range.py 150000000 150050000 [depth_k]

Survivors are merged into farm_survivors.json, the file p_hunter_eliminate.py already reads, so
afterwards just run the eliminator (it picks up the new entries). Existing survivors in that file
are kept; the swept range is added to "prior_ranges_swept". Each range also keeps a checkpoint
file range_<start>_<end>.json, so an interrupted scan resumes where it stopped.

Survivors are NOT proven prime: only Lucas-Lehmer / PRP can do that.
"""
import json
import math
import os
import sys
import time

import p_hunter_farm as farm
import revised_binV3 as core

TARGET_FILE = "farm_survivors.json"
FARM_STATE = "farm_state.json"
FLOOR_P = 136279841          # largest known Mersenne exponent
DEFAULT_START = 136279843
DEFAULT_DEPTH = 10**8        # prefilter only; the eliminator goes far deeper
CORE_SECONDS_PER_K = 5e-8    # measured p_hunter_tf rate, one core


def num(text):
    return int(float(str(text).replace(",", "").replace("_", "")))


def ask(prompt, default):
    shown = f"{default:,}" if isinstance(default, int) else default
    raw = input(f"{prompt} [{shown}]: ").strip()
    return default if raw == "" else num(raw) if isinstance(default, int) else raw


def merge_ranges(ranges):
    out = []
    for lo, hi in sorted(ranges):
        if out and lo <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return out


def range_file(start, end):
    return f"range_{start}_{end}.json"


def load_target():
    """Return (survivors by p, swept ranges, next_p) from the existing farm_survivors.json."""
    survivors, swept, next_p = {}, None, None
    if os.path.exists(TARGET_FILE):
        with open(TARGET_FILE) as f:
            src = json.load(f)
        survivors = {s["p"]: s for s in src["survivors"]}
        swept = src.get("prior_ranges_swept")
        next_p = src.get("next_p")
    if swept is None:
        swept = []
        if os.path.exists(FARM_STATE):
            with open(FARM_STATE) as f:
                fs = json.load(f)
            if fs.get("start_p") is not None:
                swept = [[fs["start_p"], fs["next_p"] - 1]]
            next_p = next_p or fs.get("next_p")
    return survivors, swept, next_p


def suggested_start(swept, next_p):
    ends = [hi + 1 for _, hi in swept] + ([next_p] if next_p else [])
    return max(ends + [DEFAULT_START])


def save_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1)
    os.replace(tmp, path)


def write_target(base, base_swept, base_next, rng):
    """Merge this range's survivors and swept span into farm_survivors.json."""
    merged = dict(base)
    for s in rng["survivors"]:
        merged.setdefault(s["p"], s)
    ranked = sorted(merged.values(), key=lambda s: (-s["bits"], -s["score"], s["p"]))
    swept = merge_ranges(base_swept + [[rng["start_p"], rng["next_p"] - 1]]) \
        if rng["next_p"] > rng["start_p"] else base_swept
    save_json(TARGET_FILE, {"next_p": base_next or rng["next_p"],
                            "prior_ranges_swept": swept, "survivors": ranked})
    return len(ranked)


def estimate(start, end, depth_k, workers):
    mid = (start + end) / 2
    primes = (end - start) / math.log(mid)
    surv = math.log(2 * mid) / math.log(2 * depth_k * mid + 1)
    secs = primes * surv * depth_k * CORE_SECONDS_PER_K / workers
    return primes, surv, secs


def fmt_time(s):
    s = int(s)
    if s >= 3600:
        return f"{s // 3600}h{s % 3600 // 60:02d}m"
    return f"{s // 60}m{s % 60:02d}s"


def scan_range(start, end, depth_k, workers, target):
    path = range_file(start, end)
    rng = None
    if os.path.exists(path):
        with open(path) as f:
            old = json.load(f)
        if old["next_p"] > end:
            print(f"{path}: already complete ({len(old['survivors'])} survivors).")
            if input("Re-run from scratch? [y/N]: ").strip().lower() != "y":
                return old
        else:
            print(f"{path}: checkpoint at p={old['next_p']:,} ({len(old['survivors'])} survivors, "
                  f"depth k<={old['depth_k']:,}).")
            if input("Resume? [Y/n]: ").strip().lower() in ("", "y", "yes"):
                rng, depth_k = old, old["depth_k"]
    if rng is None:
        rng = {"start_p": start, "end_p": end, "depth_k": depth_k, "next_p": start,
               "survivors": [], "eliminated": 0, "primes_scanned": 0}

    base, base_swept, base_next = load_target()
    chunk = 4 * workers
    found, t0 = 0, time.time()
    print(f"\nScanning p in [{rng['next_p']:,}, {end:,}]  depth k<={depth_k:,} "
          f"(~2^{farm.depth_bits(rng['next_p'], depth_k):.1f})  workers={workers}"
          + (f"  stop after {target} new survivors" if target else ""))
    pool = farm.make_pool(workers)
    try:
        while rng["next_p"] <= end and (not target or found < target):
            primes, scanned_to = farm.next_primes(rng["next_p"], chunk)
            primes = [p for p in primes if p <= end]
            scanned_to = min(scanned_to, end + 1)
            tc = time.time()
            results = list(pool.map(farm.probe_p, [(p, depth_k) for p in primes]))
            alive = sorted(p for p, hit in results if hit is None)
            killed = len(primes) - len(alive)

            if target and found + len(alive) >= target:
                alive = alive[:target - found]
                cut = alive[-1]
                kept = [p for p in primes if p <= cut]
                killed = len(kept) - len(alive)
                rng["next_p"] = cut + 1
                primes = kept
            else:
                rng["next_p"] = scanned_to

            for p in alive:
                rng["survivors"].append({"p": p, "depth_k": depth_k,
                                         "bits": round(farm.depth_bits(p, depth_k), 2),
                                         "score": core.mersenne_likeness_score(p)[0]})
            found += len(alive)
            rng["eliminated"] += killed
            rng["primes_scanned"] += len(primes)
            save_json(path, rng)
            total = write_target(base, base_swept, base_next, rng)
            done = (rng["next_p"] - start) / max(end + 1 - start, 1)
            print(f"  p<{rng['next_p']:,}  {done:6.1%}  tested={len(primes)} killed={killed} "
                  f"survived={len(alive)}  range total={len(rng['survivors'])}  "
                  f"({time.time() - tc:.1f}s, {fmt_time(time.time() - t0)} elapsed)")
    except KeyboardInterrupt:
        print("\nInterrupted; the last finished batch is saved. Run the same range again to resume.")
        pool.shutdown(wait=False, cancel_futures=True)
        return rng
    pool.shutdown()
    print(f"\nRange scan {'complete' if rng['next_p'] > end else 'stopped at target'}: "
          f"{len(rng['survivors'])} survivors from {rng['primes_scanned']:,} primes "
          f"(next unscanned p = {rng['next_p']:,}).")
    return rng


def main():
    if farm.BACKEND == "none":
        sys.exit("No Rust backend. Build p_hunter_tf: cd p_hunter_tf && maturin develop --release")
    workers = os.cpu_count() or 1
    _, swept, next_p = load_target()

    if len(sys.argv) >= 3:
        start, end = num(sys.argv[1]), num(sys.argv[2])
        depth_k = num(sys.argv[3]) if len(sys.argv) >= 4 else DEFAULT_DEPTH
        target = 0
        args_mode = True
    else:
        args_mode = False

    while True:
        if not args_mode:
            print()
            print("=" * 70)
            print("P_HUNTER RANGE SCAN")
            print("=" * 70)
            if swept:
                print("Already swept: " + ", ".join(f"{a:,}-{b:,}" for a, b in swept))
            start = ask("Start p", suggested_start(swept, next_p))
            end = ask("End p", start + 50_000)
            depth_k = ask("Prefilter depth, max k", DEFAULT_DEPTH)
            target = ask("Stop after N survivors (0 = whole range)", 0)
            workers = ask("Workers", workers)

        if start <= FLOOR_P:
            print(f"Start must be above {FLOOR_P:,} (the largest known exponent).")
        elif end < start or depth_k < 1 or workers < 1 or target < 0:
            print("Invalid range or settings.")
        elif farm.BACKEND == "p_hunter_rs" and depth_k > farm.RS_MAX_K:
            print(f"Depth above {farm.RS_MAX_K:,} is not rigorous with p_hunter_rs; build p_hunter_tf.")
        else:
            primes, surv, secs = estimate(start, end, depth_k, workers)
            print(f"\nAbout {primes:,.0f} primes; expect ~{surv:.0%} to survive k<={depth_k:,} "
                  f"(~{primes * surv:,.0f} survivors), roughly {fmt_time(secs)}"
                  + (" if the whole range is scanned." if target else "."))
            if args_mode or input("Start? [Y/n]: ").strip().lower() in ("", "y", "yes"):
                rng = scan_range(start, end, depth_k, workers, target)
                print(f"Survivors merged into {TARGET_FILE}. Next: python p_hunter_eliminate.py -> option 1.")
                _, swept, next_p = load_target()
        if args_mode:
            break
        if input("\nScan another range? [y/N]: ").strip().lower() != "y":
            break


if __name__ == "__main__":
    try:
        main()
    except (EOFError, KeyboardInterrupt):
        print()
