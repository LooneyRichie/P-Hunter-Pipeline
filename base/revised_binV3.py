"""
revised_binV3.py — fast Mersenne exponent candidate sieve.

CORRECTNESS MODEL
-----------------
Rigorous filters (zero false negatives):
  1. p ∈ {2} ∪ {odd primes}          — necessary for M_p prime
  2. no prime factor q = 2kp+1 of M_p up to max_k  — small-factor exclusion

Ranking heuristics (NOT filters, NEVER used in the candidate path):
  - mersenne_likeness_score(p)   — sorts survivors by binary shape
  - calibrated_probability(...)  — logistic fit, needs labeled data

There is NO theorem constraining the binary shape of a Mersenne exponent
beyond bit0(p) = 1. p = 2203 (7 runs) and p = 9689 (9 runs) are proven
Mersenne exponents. Any max_runs / shape gate deletes real candidates.
"""

import math
import os
import json
import time
import random
from concurrent.futures import ProcessPoolExecutor, as_completed

MR_BASES_64 = (2, 325, 9375, 28178, 450775, 9780504, 1795265022)
STATE_FILE = "p_hunter_state.json"
CALIBRATION_FILE = "scorer_calibration.json"
DEEP_STATE_FILE = "deep_verify_state.json"

SCAN_MAX_K = 1000
DEEP_MAX_K = 2000000000

KNOWN_MERSENNE_EXPONENTS = {
    2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281,
    3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243,
    110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221,
    3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457,
    32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917,
    82589933, 136279841,
}

KNOWN_COMPOSITE_DOUBLE_MERSENNE_EXPONENTS = {131071, 524287, 2147483647}

SMALL_KNOWN_P = frozenset({2, 3, 5, 7})
_WHEEL_WHITELIST = frozenset({2, 3, 5, 7, 11, 13})


# ============================================================
# STATE
# ============================================================

def load_state():
    if os.path.exists(STATE_FILE):
        try:
            with open(STATE_FILE, "r") as f:
                return json.load(f)
        except Exception:
            pass
    return None


def save_state(last_p, batch_size=100000, candidates=None):
    state = {
        "last_p": last_p,
        "batch_size": batch_size,
        "candidates": candidates or [],
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    try:
        with open(STATE_FILE, "w") as f:
            json.dump(state, f, indent=2)
    except Exception as e:
        print(f"Warning: could not save state: {e}")


def _load_deep_state():
    if not os.path.exists(DEEP_STATE_FILE):
        return None
    try:
        with open(DEEP_STATE_FILE, "r") as f:
            return json.load(f)
    except Exception:
        return None


def _save_deep_state(state):
    """Atomic-ish write: write to tmp, then rename. Prevents a crash
    mid-write from leaving a corrupt checkpoint."""
    tmp = DEEP_STATE_FILE + ".tmp"
    try:
        with open(tmp, "w") as f:
            json.dump(state, f)
        os.replace(tmp, DEEP_STATE_FILE)
    except Exception as e:
        print(f"[deep] warning: could not save checkpoint: {e}")


def _clear_deep_state():
    for f in (DEEP_STATE_FILE, DEEP_STATE_FILE + ".tmp"):
        try:
            if os.path.exists(f):
                os.remove(f)
        except Exception:
            pass


# ============================================================
# PRIMALITY (deterministic for 64-bit)
# ============================================================

def is_prime(n):
    if n < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n == sp:
            return True
        if n % sp == 0:
            return False
    d = n - 1
    r = 0
    while d % 2 == 0:
        r += 1
        d //= 2
    for base in MR_BASES_64:
        if base % n == 0:
            continue
        w = pow(base, d, n)
        if w == 1 or w == n - 1:
            continue
        for _ in range(r - 1):
            w = (w * w) % n
            if w == n - 1:
                break
        else:
            return False
    return True


# ============================================================
# GATE 0 — small-prime wheel
# ============================================================

WHEEL_MOD = 2 * 3 * 5 * 7 * 11 * 13   # 30030


def _build_wheel():
    w = bytearray(b"\x01") * WHEEL_MOD
    for m in (2, 3, 5, 7, 11, 13):
        w[0::m] = b"\x00" * (WHEEL_MOD // m)
    return bytes(w)


_WHEEL = _build_wheel()


def gate0_candidate_stream(p):
    if p < 2:
        return False
    if p <= 13:
        return p in _WHEEL_WHITELIST
    return _WHEEL[p % WHEEL_MOD] == 1


# ============================================================
# GATE A1 — rigorous binary signature
# ============================================================

def gateA1_binary_signature(p):
    if p == 2:
        return True, {"bit0": 1, "special": "p=2"}
    if p < 2 or (p & 1) == 0:
        return False, {"reason": "not_odd"}
    return True, {"bit0": 1}


# ============================================================
# GATE B2 — small-factor probe on M_p
# ============================================================

_Q_SCREEN = (5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
_SMALL_MR_BASES = (2, 3, 5)
_SMALL_MR_LIMIT = 25_326_001
_Q_RESIDUES = (1, 7, 17, 23)


def _is_prime_small(n):
    if n < 2:
        return False
    d = n - 1
    r = 0
    while (d & 1) == 0:
        r += 1
        d >>= 1
    for base in _SMALL_MR_BASES:
        if base % n == 0:
            continue
        w = pow(base, d, n)
        if w == 1 or w == n - 1:
            continue
        for _ in range(r - 1):
            w = (w * w) % n
            if w == n - 1:
                break
        else:
            return False
    return True


def gateB2_factor_probe(p, max_k=SCAN_MAX_K):
    if p in SMALL_KNOWN_P:
        return True, {"small_p_known": True}

    mp = (1 << p) - 1
    q_max = 2 * max_k * p + 1
    prime_check = _is_prime_small if q_max < _SMALL_MR_LIMIT else is_prime

    for k in range(1, max_k + 1):
        q = 2 * k * p + 1
        if q >= mp:            # <-- NEW: stop before q = mp
            break
        if q % 24 not in _Q_RESIDUES:
            continue

        survived = True
        for r in _Q_SCREEN:
            if r >= q:
                break
            if q % r == 0:
                survived = False
                break
        if not survived:
            continue

        if not prime_check(q):
            continue
        if pow(2, p, q) == 1:
            return False, {"factor": q, "k": k}

    return True, {"max_k": max_k}

def gateB3_deep_factor_probe(p, max_k=DEEP_MAX_K):
    return gateB2_factor_probe(p, max_k=max_k)


# ============================================================
# GATE C — perfect-number last-digit check
# ============================================================

def gateC2_last_digit(p):
    a = pow(2, p - 1, 10)
    b = pow(2, p, 10)
    d = (a * (b - 1)) % 10
    return d in (6, 8), {"N_mod_10": d}


# ============================================================
# PIPELINE
# ============================================================

def evaluate_candidate(p):
    if p in KNOWN_COMPOSITE_DOUBLE_MERSENNE_EXPONENTS:
        return False
    if not gate0_candidate_stream(p):
        return False
    ok, _ = gateA1_binary_signature(p)
    if not ok:
        return False
    if p in SMALL_KNOWN_P:
        return True
    if not is_prime(p):
        return False
    ok, _ = gateB2_factor_probe(p, max_k=SCAN_MAX_K)
    if not ok:
        return False
    ok, _ = gateC2_last_digit(p)
    return ok


# ============================================================
# RANKING-ONLY heuristics
# ============================================================

def _runs_of(p):
    b = bin(p)[2:]
    out = []
    cur, ln = b[0], 1
    for ch in b[1:]:
        if ch == cur:
            ln += 1
        else:
            out.append((cur, ln))
            cur, ln = ch, 1
    out.append((cur, ln))
    return out


def mersenne_likeness_score(p):
    if p < 2:
        return 0, {}
    runs = _runs_of(p)
    one_runs = [l for c, l in runs if c == '1']
    zero_runs = [l for c, l in runs if c == '0']
    if not one_runs:
        return 0, {}
    breakdown = {}
    if one_runs[0] % 2 == 1:
        breakdown["first_1run_odd"] = 2
    if one_runs[-1] % 2 == 1:
        breakdown["last_1run_odd"] = 3
    longest_one = max(one_runs)
    longest_zero = max(zero_runs) if zero_runs else 0
    if longest_zero >= longest_one:
        breakdown["long_zero_dominates"] = 3
    if sum(1 for l in one_runs if l > 2) <= 1:
        breakdown["single_long_1run"] = 2
    return sum(breakdown.values()), breakdown


def score_survivors(survivors):
    scored = []
    for p in survivors:
        s, b = mersenne_likeness_score(p)
        scored.append((s, p, b))
    scored.sort(key=lambda t: (-t[0], t[1]))
    return scored


def validate_likeness_features(start_p=3, end_p=200_000, negative_sample=5000):
    pos = [p for p in KNOWN_MERSENNE_EXPONENTS if start_p <= p <= end_p]
    pool = [p for p in range(start_p, end_p + 1)
            if p not in KNOWN_MERSENNE_EXPONENTS and is_prime(p)]
    if not pool:
        print("[validate] no negatives in range.")
        return
    neg = random.sample(pool, min(negative_sample, len(pool)))
    keys = ("first_1run_odd", "last_1run_odd",
            "long_zero_dominates", "single_long_1run")

    def rate(ps):
        c = {k: 0 for k in keys}
        for p in ps:
            _, b = mersenne_likeness_score(p)
            for k in b:
                if k in c:
                    c[k] += 1
        n = max(len(ps), 1)
        return {k: v / n for k, v in c.items()}

    fp, fn = rate(pos), rate(neg)
    print(f"{'feature':24s}  {'pos':>7s}  {'neg':>7s}  {'delta':>8s}")
    for k in keys:
        d = fp[k] - fn[k]
        flag = "" if abs(d) >= 0.05 else "  <-- noise, drop"
        print(f"{k:24s}  {fp[k]:>7.3f}  {fn[k]:>7.3f}  {d:>+8.3f}{flag}")


# ============================================================
# CALIBRATION
# ============================================================

def load_calibration():
    if not os.path.exists(CALIBRATION_FILE):
        return None
    try:
        with open(CALIBRATION_FILE, "r") as f:
            return json.load(f)
    except Exception:
        return None


def calibrated_probability(features, weights):
    logit = weights.get("intercept", 0.0)
    for k, v in features.items():
        logit += weights.get(k, 0.0) * v
    return 1.0 / (1.0 + math.exp(-logit))


def heuristic_density_probability(p):
    if p < 3:
        return 1.0
    return 1.0 / (p * math.log(2))


def build_calibration_set(start_p=3, end_p=200_000, negative_sample=5000):
    pos = [p for p in KNOWN_MERSENNE_EXPONENTS if start_p <= p <= end_p]
    pool = [p for p in range(start_p, end_p + 1)
            if p not in KNOWN_MERSENNE_EXPONENTS and is_prime(p)]
    neg = random.sample(pool, min(negative_sample, len(pool)))

    def featurize(p):
        return {
            "passed_gate0": int(gate0_candidate_stream(p)),
            "passed_A1":    int(gateA1_binary_signature(p)[0]),
            "passed_prime": int(is_prime(p)),
            "passed_B2":    int(gateB2_factor_probe(p, max_k=SCAN_MAX_K)[0]),
            "passed_C2":    int(gateC2_last_digit(p)[0]),
        }

    X = [featurize(p) for p in pos + neg]
    y = [1] * len(pos) + [0] * len(neg)
    return X, y


def fit_and_save_calibration(start_p=3, end_p=200_000, negative_sample=5000):
    try:
        from sklearn.linear_model import LogisticRegression
    except ImportError:
        print("[calib] needs scikit-learn: pip install scikit-learn")
        return None
    X, y = build_calibration_set(start_p, end_p, negative_sample)
    if not X or sum(y) == 0:
        print("[calib] empty or single-class set; widen range.")
        return None
    keys = sorted({k for row in X for k in row})
    Xm = [[row.get(k, 0) for k in keys] for row in X]
    clf = LogisticRegression(max_iter=1000, class_weight="balanced")
    clf.fit(Xm, y)
    weights = {"intercept": float(clf.intercept_[0])}
    weights.update({k: float(c) for k, c in zip(keys, clf.coef_[0])})
    with open(CALIBRATION_FILE, "w") as f:
        json.dump(weights, f, indent=2)
    print(f"[calib] saved {CALIBRATION_FILE}")
    return weights


# ============================================================
# ANALYSIS / REPORT
# ============================================================

def analyze_p(p, max_k=SCAN_MAX_K):
    print()
    print("=" * 70)
    print(f"ANALYSIS FOR p = {p}")
    print("=" * 70)
    print(f"binary: {bin(p)[2:]}")

    g0 = gate0_candidate_stream(p)
    g1, g1info = gateA1_binary_signature(p)
    pr = is_prime(p) if (g0 and g1) else False
    gB, gBinfo = (gateB2_factor_probe(p, max_k=max_k)
                  if (g0 and g1 and pr) else (False, {"skipped": True}))
    gC, gCinfo = (gateC2_last_digit(p)
                  if (g0 and g1 and pr and gB) else (False, {"skipped": True}))

    print(f"\nGate 0  wheel          : {'PASS' if g0 else 'FAIL'}")
    print(f"Gate A1 bit0           : {'PASS' if g1 else 'FAIL'} | {g1info}")
    print(f"Gate 2  is_prime(p)    : {'PASS' if pr else 'FAIL'}")
    print(f"Gate B2 2kp+1 probe    : {'PASS' if gB else 'FAIL'} | {gBinfo}")
    print(f"Gate C2 last digit     : {'PASS' if gC else 'FAIL'} | {gCinfo}")

    print("\n--- RANKING (not filters, not probabilities) ---")
    score, brk = mersenne_likeness_score(p)
    print(f"Likeness score: {score}/10  features={list(brk.keys())}")

    print("\n--- HEURISTIC PROBABILITY ---")
    print(f"Baseline 1/(p·ln2): {heuristic_density_probability(p):.3e}")
    calib = load_calibration()
    if calib is not None:
        feats = {
            "passed_gate0": int(g0),
            "passed_A1":    int(g1),
            "passed_prime": int(pr),
            "passed_B2":    int(gB),
            "passed_C2":    int(gC),
        }
        print(f"Calibrated P(M_p prime | features): "
              f"{calibrated_probability(feats, calib):.4%}")
    else:
        print("Calibrated probability: not fitted. "
              "Run fit_and_save_calibration() once.")

    verdict = g0 and g1 and pr and gB and gC
    print("\nVERDICT: " + ("SURVIVOR (no small factor found)" if verdict
                            else "ELIMINATED"))
    if verdict:
        print("  Not a primality proof. Only Lucas-Lehmer proves M_p prime.")
    print("=" * 70)


# ============================================================
# SEARCH DRIVERS
# ============================================================

def run_range_search(start_p=0, end_p=10000, verbose=True):
    print()
    print("=" * 70)
    print(f"RANGE SEARCH: p ∈ [{start_p}, {end_p}]")
    print("=" * 70)
    t0 = time.time()
    survivors = [p for p in range(start_p, end_p + 1) if evaluate_candidate(p)]
    elapsed = time.time() - t0
    total = end_p - start_p + 1

    print(f"Range size : {total:,}")
    print(f"Survivors  : {len(survivors)}")
    print(f"Time       : {elapsed:.2f}s")
    print(f"Throughput : {total/elapsed:,.0f} nums/sec")
    if survivors and verbose:
        print("\nTop survivors by likeness score:")
        for s, p, b in score_survivors(survivors)[:50]:
            tag = " [KNOWN]" if p in KNOWN_MERSENNE_EXPONENTS else ""
            print(f"  score={s:>2}/10  p={p:<12}  bin={bin(p)[2:]}{tag}")
    print("=" * 70)
    return survivors


def run_benchmark(start_p=0, end_p=10000):
    known = sorted(p for p in KNOWN_MERSENNE_EXPONENTS if start_p <= p <= end_p)
    print()
    print("=" * 70)
    print(f"BENCHMARK p ∈ [{start_p}, {end_p}]")
    print("=" * 70)
    print(f"Known Mersenne exponents in range: {known}")

    t0 = time.time()
    survivors = [p for p in range(start_p, end_p + 1) if evaluate_candidate(p)]
    elapsed = time.time() - t0

    captured = sorted(p for p in survivors if p in KNOWN_MERSENNE_EXPONENTS)
    unmatched = sorted(p for p in survivors if p not in KNOWN_MERSENNE_EXPONENTS)
    missed = sorted(p for p in known if p not in survivors)

    print(f"\nRange size             : {end_p - start_p + 1:,}")
    print(f"Survivors              : {len(survivors)}")
    print(f"Known captured         : {len(captured)}/{len(known)}")
    if missed:
        print(f"  MISSED (should be 0) : {missed}")
    print(f"Unmatched survivors    : {len(unmatched)}")
    print(f"  (survivors with no known factor ≤ 2·{SCAN_MAX_K}·p+1 — not false positives)")
    print(f"Time                   : {elapsed:.2f}s")
    print("=" * 70)
    return survivors


def run_endless_search(start_p=0, batch_size=100000, candidates_list=None):
    current_p = start_p
    candidates = candidates_list or []
    print()
    print("=" * 70)
    print(f"ENDLESS SEARCH from p = {current_p:,}  |  batch = {batch_size:,}")
    print("Ctrl+C to checkpoint and exit.")
    print("=" * 70)
    try:
        while True:
            batch_end = current_p + batch_size - 1
            print(f"--> [{current_p:,} .. {batch_end:,}]")
            found = False
            for p in range(current_p, batch_end + 1):
                if evaluate_candidate(p):
                    s, b = mersenne_likeness_score(p)
                    candidates.append(p)
                    save_state(p, candidates=candidates, batch_size=batch_size)
                    print()
                    print("=" * 70)
                    print(f"CANDIDATE p = {p:,}  score={s}/10  feats={list(b.keys())}")
                    print(f"binary: {bin(p)[2:]}")
                    print("=" * 70)
                    ch = input("[C]ontinue / [A]nalyze / [Q]uit [C]: ").strip().lower()
                    if ch == "a":
                        analyze_p(p)
                        input("Enter to resume...")
                    elif ch == "q":
                        return candidates
                    current_p = p + 1
                    found = True
                    break
            if not found:
                current_p = batch_end + 1
                save_state(current_p - 1, candidates=candidates,
                           batch_size=batch_size)
                print(f"   batch complete, saved at {current_p - 1:,}")
    except KeyboardInterrupt:
        save_state(current_p, candidates=candidates, batch_size=batch_size)
        print(f"\nCheckpoint saved at p = {current_p:,}")
    return candidates


def run_continue_search():
    st = load_state()
    if not st:
        print("[!] no checkpoint. Starting fresh.")
        return run_endless_search(0)
    print()
    print("=" * 70)
    print(f"RESUME  last_p={st.get('last_p', 0):,}  "
          f"candidates={len(st.get('candidates', []))}")
    print(f"saved at {st.get('timestamp', '?')}")
    print("=" * 70)
    nxt = st.get("last_p", 0) + 1
    if input(f"Resume from {nxt:,}? [Y/n]: ").strip().lower() in ("", "y", "yes"):
        return run_endless_search(nxt,
                                  batch_size=st.get("batch_size", 100000),
                                  candidates_list=st.get("candidates", []))
    return st.get("candidates", [])


def _scan_chunk(args):
    a, b = args
    return [p for p in range(a, b + 1) if evaluate_candidate(p)]


def _chunks(a, b, size):
    while a <= b:
        e = min(a + size - 1, b)
        yield (a, e)
        a = e + 1


def run_parallel_range_search(start_p=0, end_p=1000000, workers=None, verbose=True):
    workers = workers or os.cpu_count() or 4
    total = end_p - start_p + 1
    chunk = max(50_000, total // (workers * 4))
    print()
    print("=" * 70)
    print(f"PARALLEL p ∈ [{start_p:,}, {end_p:,}]  workers={workers}  chunk={chunk:,}")
    print("=" * 70)
    t0 = time.time()
    survivors = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_scan_chunk, c) for c in _chunks(start_p, end_p, chunk)]
        for f in as_completed(futs):
            survivors.extend(f.result())
    survivors.sort()
    elapsed = time.time() - t0
    print(f"Range      : {total:,}")
    print(f"Survivors  : {len(survivors)}")
    print(f"Time       : {elapsed:.2f}s")
    print(f"Throughput : {total/elapsed:,.0f} nums/sec")
    if verbose and survivors:
        print("\nTop survivors:")
        for s, p, b in score_survivors(survivors)[:50]:
            tag = " [KNOWN]" if p in KNOWN_MERSENNE_EXPONENTS else ""
            print(f"  score={s:>2}/10  p={p:<12}  bin={bin(p)[2:]}{tag}")
    print("=" * 70)
    return survivors


# ============================================================
# DEEP VERIFY WITH CHECKPOINTING
# ============================================================

def run_deep_verify(start_p, end_p, survivors, max_k=DEEP_MAX_K,
                    verbose=True, autosave_every=10, autosave_secs=5.0):
    """
    Deep factor probe over survivors, with checkpoint/resume.

    State file: DEEP_STATE_FILE. Written every `autosave_every` survivors
    OR every `autosave_secs` seconds (whichever first), and on Ctrl+C.

    On resume, the survivor list from the checkpoint is used — not the
    one passed in — so the resumed run is bit-identical to what it would
    have been without the interruption.

    Ctrl+C saves and returns partial results without raising.
    """
    existing = _load_deep_state()
    resuming = False

    if existing is not None:
        if (existing.get("start_p") == start_p and
            existing.get("end_p") == end_p and
            existing.get("max_k") == max_k):
            total_done = existing.get("processed_idx", 0)
            total_surv = len(existing.get("survivors", []))
            print()
            print(f"[deep] checkpoint found: {total_done}/{total_surv} survivors "
                  f"processed, saved {existing.get('timestamp', '?')}")
            ans = input("Resume from checkpoint? [Y/n] ").strip().lower()
            if ans in ("", "y", "yes"):
                resuming = True
                survivors = existing["survivors"]
            else:
                _clear_deep_state()
        else:
            print(f"[deep] checkpoint exists for a different run "
                  f"({existing.get('start_p')}..{existing.get('end_p')}, "
                  f"max_k={existing.get('max_k')}). Leaving it alone.")
            print(f"[deep] Starting fresh; existing checkpoint will be "
                  f"overwritten on first save.")

    if resuming:
        state = existing
    else:
        scored = score_survivors(survivors)
        survivors_sorted = [p for _, p, _ in scored]
        state = {
            "start_p": start_p,
            "end_p": end_p,
            "max_k": max_k,
            "survivors": survivors_sorted,
            "processed_idx": 0,
            "known": [],
            "eliminated": [],
            "unverified": [],
            "complete": False,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        _save_deep_state(state)

    surv_list = state["survivors"]
    total = len(surv_list)
    idx = state["processed_idx"]
    known = list(state["known"])
    eliminated = [tuple(x) for x in state["eliminated"]]
    unverified = list(state["unverified"])

    print()
    print("=" * 70)
    print(f"DEEP VERIFY  survivors={total}  max_k={max_k:,}  "
          f"start_idx={idx}")
    print(f"Checkpoint : {DEEP_STATE_FILE}  "
          f"(every {autosave_every} survivors or {autosave_secs:.0f}s)")
    print("Ctrl+C saves and exits cleanly.")
    print("=" * 70)

    t0 = time.time()
    last_save = t0
    since_save = 0

    def _flush():
        state["processed_idx"] = idx
        state["known"] = known
        state["eliminated"] = [list(x) for x in eliminated]
        state["unverified"] = unverified
        state["timestamp"] = time.strftime("%Y-%m-%d %H:%M:%S")
        _save_deep_state(state)

    try:
        while idx < total:
            p = surv_list[idx]
            if p in KNOWN_MERSENNE_EXPONENTS:
                known.append(p)
            else:
                ok, info = gateB3_deep_factor_probe(p, max_k=max_k)
                if ok:
                    unverified.append(p)
                else:
                    eliminated.append((p, info["factor"]))
                    if verbose:
                        print(f"[-] p={p}  factor q={info['factor']}")
            idx += 1
            since_save += 1

            now = time.time()
            if since_save >= autosave_every or now - last_save >= autosave_secs:
                _flush()
                last_save = now
                since_save = 0
                if verbose:
                    pct = idx / total
                    elapsed = now - t0
                    print(f"    ... {idx}/{total} ({pct:.1%})  "
                          f"elapsed={elapsed:.0f}s  "
                          f"elim={len(eliminated)}  "
                          f"unver={len(unverified)}")

        # Complete
        state["complete"] = True
        _flush()

    except KeyboardInterrupt:
        _flush()
        print()
        print(f"[deep] interrupted at {idx}/{total}. "
              f"Checkpoint saved to {DEEP_STATE_FILE}.")
        print(f"[deep] Run menu 6 again with the same range to resume.")
        return {"known": known, "eliminated": eliminated,
                "unverified": unverified,
                "interrupted": True, "idx": idx, "total": total}

    elapsed = time.time() - t0
    print()
    print(f"Time         : {elapsed:.2f}s")
    print(f"Known        : {known}")
    print(f"Eliminated   : {len(eliminated)}")
    print(f"Unverified   : {len(unverified)}")
    print(f"  (unverified means no factor ≤ 2·{max_k}·p+1 — NOT proven prime)")
    if unverified:
        shown = unverified[:20]
        print(f"  Unverified list (first 20): {shown}")
    print("=" * 70)
    print(f"[deep] complete. Checkpoint retained at {DEEP_STATE_FILE}.")
    print(f"[deep] Delete that file, or start a new range, to clear.")

    return {"known": known, "eliminated": eliminated,
            "unverified": unverified,
            "interrupted": False, "idx": total, "total": total}


def run_range_search_and_deep_verify(start_p, end_p, max_k=DEEP_MAX_K):
    """
    Combined menu-6 flow with deep-verify checkpointing.

    Detects a matching checkpoint and offers to skip range search and
    resume deep verify directly. If no checkpoint matches, runs the
    full flow: range search, then deep verify.
    """
    existing = _load_deep_state()
    if (existing is not None and
        existing.get("start_p") == start_p and
        existing.get("end_p") == end_p and
        existing.get("max_k") == max_k and
        not existing.get("complete", False)):
        total_surv = len(existing.get("survivors", []))
        done = existing.get("processed_idx", 0)
        print()
        print(f"[deep] matching checkpoint for p ∈ [{start_p}, {end_p}], "
              f"max_k={max_k}: {done}/{total_surv} done")
        ans = input("Resume deep verify without re-running range search? "
                    "[Y/n] ").strip().lower()
        if ans in ("", "y", "yes"):
            survivors = existing["survivors"]
            return run_deep_verify(start_p, end_p, survivors, max_k=max_k)

    survivors = run_range_search(start_p, end_p, verbose=False)
    return run_deep_verify(start_p, end_p, survivors, max_k=max_k)


# ============================================================
# MENU
# ============================================================

def main_menu():
    while True:
        calib = "fitted" if os.path.exists(CALIBRATION_FILE) else "not fitted"
        deep_pending = _load_deep_state()
        if deep_pending and not deep_pending.get("complete", False):
            done = deep_pending.get("processed_idx", 0)
            tot = len(deep_pending.get("survivors", []))
            deep_status = (f"paused at {done}/{tot} "
                           f"({deep_pending.get('start_p')}.."
                           f"{deep_pending.get('end_p')})")
        else:
            deep_status = "none"

        print()
        print("=" * 70)
        print("P_HUNTER")
        print("=" * 70)
        print(f"Pipeline : wheel -> bit0 -> is_prime -> 2kp+1 probe -> mod 10")
        print(f"Calib    : {calib}")
        print(f"Deep     : {deep_status}")
        print("-" * 70)
        print("1. Range search")
        print("2. Endless search")
        print("3. Resume last session")
        print("4. Benchmark vs known exponents")
        print("5. Analyze single p")
        print("6. Range search + deep verify  (checkpointed)")
        print("7. Parallel range search")
        print("8. Validate likeness features")
        print("9. Fit calibration (needs scikit-learn)")
        print("D. Discard deep-verify checkpoint")
        print("0. Exit")
        print("=" * 70)
        c = input("Choice [0-9,D]: ").strip()

        try:
            if c == "1":
                a = int(input("start [0]: ") or "0")
                b = int(input("end [10000]: ") or "10000")
                run_range_search(a, b)
            elif c == "2":
                a = int(input("start [0]: ") or "0")
                bs = int(input("batch [100000]: ") or "100000")
                run_endless_search(a, bs)
            elif c == "3":
                run_continue_search()
            elif c == "4":
                a = int(input("start [0]: ") or "0")
                b = int(input("end [10000]: ") or "10000")
                run_benchmark(a, b)
            elif c == "5":
                p = int(input("p: "))
                mk = int(input(f"max_k [{SCAN_MAX_K}]: ") or str(SCAN_MAX_K))
                analyze_p(p, max_k=mk)
            elif c == "6":
                a = int(input("start [0]: ") or "0")
                b = int(input("end [10000]: ") or "10000")
                mk = int(input(f"deep max_k [{DEEP_MAX_K}]: ") or str(DEEP_MAX_K))
                run_range_search_and_deep_verify(a, b, max_k=mk)
            elif c == "7":
                a = int(input("start [0]: ") or "0")
                b = int(input("end [1000000]: ") or "1000000")
                w = input("workers [all]: ").strip()
                run_parallel_range_search(a, b, workers=int(w) if w else None)
            elif c == "8":
                a = int(input("start [3]: ") or "3")
                b = int(input("end [200000]: ") or "200000")
                validate_likeness_features(a, b)
            elif c == "9":
                a = int(input("start [3]: ") or "3")
                b = int(input("end [200000]: ") or "200000")
                fit_and_save_calibration(a, b)
            elif c.lower() == "d":
                if os.path.exists(DEEP_STATE_FILE):
                    if input("Delete deep-verify checkpoint? [y/N]: ").strip().lower() == "y":
                        _clear_deep_state()
                        print("[deep] checkpoint deleted.")
                else:
                    print("[deep] no checkpoint to delete.")
            elif c == "0":
                print("bye.")
                break
            else:
                print("invalid.")
        except ValueError:
            print("invalid input.")


if __name__ == "__main__":
    main_menu()