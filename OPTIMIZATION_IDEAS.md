# Optimization Ideas for `p_hunter_tf` / `p_hunter_rs`

A design document, not a code change. Every idea has: **why**, **estimated gain**, **steps**, **validation**, **risk**.
Follow the rules in [CONTEXT.md](CONTEXT.md): copy a file before risky edits, benchmark before and after, and **zero missed known factors is non-negotiable**.

---

## 0. Baseline: where the time goes

Hardware probed: Intel N200, 4 cores, AVX2/BMI2/ADX, **no AVX-512**, L1d 32 KiB/core, L2 2 MiB, L3 6 MiB, no NVIDIA GPU (Intel UHD iGPU only). Rust 1.89, no `.cargo/config.toml` (so no `target-cpu=native`).

Single-core `tf_window` at p=136279843, 5e7-wide bands:

| k band | q size | ns per k |
|---|---|---|
| ≤5e7 | ~2^53 | 11.3 (a factor was found early) |
| ~1e9 | ~2^58 | 41.9 |
| ~5e10 | <2^64 | 43.2 |
| ~1e11 | >2^64 | **201.9** |
| ~1e13 | ~2^71 | **210.8** |

Findings:
1. Crossing 2^64 makes each k about **4.7× slower**. `mulmod_big` uses 24-bit chunks and many u128 `%` operations (software division). This is the deep-TF bottleneck.
2. After sieving to 20M, ~3.3% of k survive; the `q mod 8 ∈ {1,7}` filter leaves ~1.7%. Each survivor costs ~27 squarings (p ≈ 2^27) plus a modular reduction each.
3. Only ~66% of sieve survivors are prime, so sieving further buys little.
4. Block size 2^23 bytes exceeds L2, so the sieve is memory-bound.
5. `inv_mod` (extended Euclid) is run for every one of ~1.27M primes in every call.
6. `p_hunter_rs` is a weaker copy: 64-bit q only (silently skips q>2^64), allocates `alive[k_max+1]`, clones the prime Vec every call, SIEVE_CAP 50M.

---

## Tier A: Kernel arithmetic (biggest wins)

### A1. Montgomery multiplication for q < 2^64 (est. 5–10×)
**Why:** `x*x % q` on u128 compiles to a `__umodti3` call (~30–100 cycles). Montgomery replaces it with 2–3 multiplications.
**Steps**
1. Copy `p_hunter_tf/src/lib.rs` to `lib.rs.bak`.
2. Precompute per q: `qinv = -q^{-1} mod 2^64` (Newton: `x = x*(2 - q*x)` 5–6 times), `r2 = 2^128 mod q`.
3. Implement `mont_mul(a,b,q,qinv)`: `t = a*b` (u128); `m = (t as u64).wrapping_mul(qinv)`; `u = (t + m*q) >> 64`; conditional subtract. Use u128 with care for overflow when q ≥ 2^63 (use the carry-aware variant).
4. Convert 1 into Montgomery form (`r mod q`), do left-to-right square-and-double. Doubling by 2 is just `x+x` with conditional subtract — no multiplication.
5. Test `result == 1` in Montgomery form (compare against `R mod q`).
6. Use only when q < 2^64; leave the old path as fallback.
**Validation:** `crosscheck.py` against brute-force `pow(2,p,q)`; random q (odd) tests; all known Mersenne primes (2^p−1 for p in known list) must never yield a factor; identical factor lists on 100 random bands.
**Risk:** low.

### A2. Multi-limb Montgomery for 64 < q < 2^128 (est. 8–20× for the deep regime)
**Why:** This removes the 4.7× cliff. Needed for any TF past 2^64 (and mandatory for the ≥2^77 plan in C2).
**Steps**
1. Represent residues as `(lo: u64, hi: u64)`; q < 2^96 keeps carries simple (MAX_Q_BITS=100 is already close).
2. Implement 2-limb CIOS Montgomery: 4 `mul` (64×64→128) for the product, 2 for reduction `m`, and 2×2 for `m*q`. Use `u64::carrying_mul`-style helpers or `core::arch::x86_64::_mulx_u64` + `_addcarryx_u64` (BMI2/ADX present).
3. Keep MAX_Q_BITS check; raise to 124 only after tests with u128 headroom.
4. Square-and-double again (doubling = shift + conditional subtract).
**Validation:** compare with Python big-int `pow(2,p,q)` on 10^6 random q in [2^64, 2^100]; known 70–90 bit factors from the GIMPS factor list must be found.
**Risk:** medium (carry bugs). Mitigate with property tests.

### A3. Skip the first squarings (est. 15–20% of modexp)
**Why:** For the top bits of p, 2^e for e < 64 is a plain shift: `2^(p >> j)` where `p >> j < 64` (fits for q > 2^63, otherwise use `< log2(q)`).
**Steps**
1. Find the largest j with `(p >> j) < bits(q)−1`.
2. Start the exponent chain from `r = 1 << (p >> j)` (already reduced since < q), process the remaining j bits.
**Validation:** same as A1. **Risk:** trivial.

### A4. Interleave candidates for instruction-level parallelism (est. 1.5–2×)
**Why:** One modexp is a serial dependency chain (each square needs the last). The CPU's multiplier sits idle waiting on latency.
**Steps**
1. Collect survivors into a batch of 4–8 candidate q's.
2. Run their square-and-double loops in lock-step (same p ⇒ same exponent bit pattern, so control flow is identical).
3. Write branch-free (`cmov`) conditional subtraction.
4. Optional later: AVX2 `vpmuludq` 4-lane 32-bit-limb Montgomery (no IFMA on N200, so scalar `mulx` is probably as fast; benchmark before committing).
**Validation:** identical results to the scalar version. **Risk:** low–medium.

### A5. Compiler settings (est. 5–20%, nearly free)
**Steps**
1. Add `.cargo/config.toml` with `rustflags = ["-C","target-cpu=native"]` (or `x86-64-v3` for portability). `lto` and `codegen-units = 1` are already set in `Cargo.toml`; optionally add `panic = "abort"`.
2. Optional: PGO (`cargo pgo`) using a representative band.
3. Rebuild with `maturin develop --release`; confirm no regressions with `crosscheck.py`.
**Risk:** a native build won't run on other CPUs.

### A6. Early-abort and cheaper test shape
Since the test is `2^p ≡ 1 (mod q)`, the result is checked only at the end; there is no early exit. Keep as is — but ensure the final comparison costs nothing in Montgomery form (compare to `R mod q`).

---

## Tier B: Sieve and wheel

### B1. k-class wheel (mfaktc style) (est. 3–5× fewer k to sieve)
**Why:** q = 2kp+1 must be ≡ ±1 (mod 8) and not divisible by 3, 5, 7, 11. That fixes which residues of k mod 4·3·5·7·11 = 4620 are possible (≈960 of 4620 classes ≈ 21%). Sieve only those classes.
**Steps**
1. For a given p, compute the allowed class list `C = {c in [0,4620): q(c) mod 8 ∈ {1,7} and q(c) mod m != 0 for m in {3,5,7,11}}`.
2. For each class, iterate `k = c + 4620·t` and sieve in `t` space (the sieve array shrinks ~4.8×).
3. Sieve primes start at 13. Each prime's `next` index is computed per class (cheap with the per-prime inverse of `4620·2p`).
4. Optionally extend to 13 (60060 → ~11,520 classes).
**Validation:** surviving k set must equal the existing sieve's survivors (set equality test on 1e7-wide windows, then on random bands).
**Risk:** medium. It's the most complex change, but a proven design.

### B2. Bitset + cache-sized segments (est. 1.5–3× sieve speed)
**Steps**
1. Replace the 8 MB byte array with a 32–64 KB (L1) or 256 KB (L2) segment.
2. Use a bitset (8× less memory traffic), or keep bytes at L1-sized segments for simplicity.
3. Keep a per-prime "next offset" array across segments so each segment resumes where the previous stopped.
**Risk:** low.

### B3. Bucket sieve for large primes (est. up to 2× on the sieve)
Primes larger than the segment size hit a segment at most once. Put them into per-segment buckets instead of scanning them every segment.
**Steps:** pre-pass pushes `(prime, offset)` into bucket[segment]; consume the bucket when processing the segment; re-bucket the next hit.

### B4. Re-tune the sieve depth after A1/A2
When modexp gets 5–20× cheaper, the optimal sieve limit drops (the sieve costs more per removed candidate than a modexp does). **Steps:** benchmark SIEVE_CAP ∈ {1M, 2M, 5M, 10M, 20M} on 3 bands (53, 60, 70-bit); pick the minimum ns/k per band; make it bit-size dependent.

### B5. Persistent per-p sieve session (est. removes ~5–15% fixed overhead per call)
**Why:** Every `tf_window` call recomputes ~1.27M extended-Euclid inversions. When bands for one p are called in sequence, those are redundant.
**Steps:** cache `(p, [inv(2p mod prime)])` in a `Mutex<HashMap>` keyed by p; compute `next[i] = -inv·... mod prime` from it. Or batch inversion (Montgomery trick: 1 inversion + 3(n−1) multiplications).

### B6. Parallelism hygiene
1. Release the GIL in the PyO3 function (`py.allow_threads`).
2. Share the prime table across threads (it's `PRIME_CACHE` already; verify there's no per-thread rebuild).
3. Use `workers = 4` (physical cores) and split bands, not primes.
4. Rayon inside the crate across classes (after B1) is a cleaner alternative to Python multiprocessing.

### B7. Retire `p_hunter_rs`
Port the sieve improvements once into `p_hunter_tf`, then have the Python wrappers use it exclusively. It is 64-bit only and silently skips q>2^64, which is a hazard.

---

## Tier C: Scheduling, economics, and "start at a certain depth"

### C1. Per-p bit-level ledger (your "list/amount idea")
**Why:** Today progress is stored as k-ranges. GIMPS tracks TF as a **bit level** b ("no factor below 2^b"). Converting makes every record comparable and mergeable with external data.
**Steps**
1. Define `k_start(b) = ceil((2^b − 1) / (2p))` and `k_end(b)`. The band for level b is `[k_start(b), k_start(b+1))`.
2. Add to `eliminate_state.json` per p: `tf_bits: 65`, and `skip_below_bits: N` meaning "start at 2^N".
3. A new level costs `2^(b−1)` k's in total (doubling per bit) — scheduling becomes "advance p to the next bit".
4. Start k for a fresh p is `k_start(skip_below_bits)` rather than 1: your "start all huge exponents at a certain depth" idea.
**Validation:** `k_start(b)` must satisfy `2·k·p+1 ≥ 2^b` and be minimal; unit-test boundary values.

### C2. Import the known TF level from PrimeNet (key finding) (est. skips weeks of redundant work)
**Evidence:** mersenne.org reports p=136282409 as *"Composite: proven not-prime with PRP test. No factors below 2^77. TF target 2^77."* Exponents in this range are typically already TF'd to ~2^77 by GIMPS GPUs, and many are PRP-tested. TF to 2^60–2^65 on this CPU is therefore very likely **redundant for the exponents you hold**.
**Steps**
1. For each survivor p, query the PrimeNet exponent report (`https://www.mersenne.org/report_exponent/?exp_lo=P&full=1`) — fetchable per exponent. (mersenne.ca per-exponent pages returned HTTP 401 for scripted access during research.)
2. Bulk alternative: https://www.mersenne.ca/export/ and PrimeNet's `/primenet/` data files list known factors and status. Download once, build a local `dict p → (tf_bits, composite?, known_factors)`.
3. Drop every p that is already composite (factor known or PRP/LL-composite) — it cannot be a Mersenne prime.
4. For the rest, set `skip_below_bits = max(known TF bits)` (C1).
5. Be polite: cache responses, ≤1 request/second, identify with a User-Agent.
**Result:** the worklist shrinks to exponents that are *unresolved*, and every CPU-second goes to bit levels nobody has covered.

### C3. Reverse the search: worklist of exponents that still need work
Instead of scanning the neighbourhood of p≈136.28M, ask: *which prime exponents have no PRP/LL result yet and TF below 2^77?* Use the same bulk data, sorted by lowest TF level / lack of tests. Those are the only places where new TF or a primality test is useful.

### C4. Breadth-first by bit level
**Why:** Per-candidate success is ≈ 1/(2k) and kill probability per bit level ≈ 1/b (Mertens). The next bit costs 2× the previous but kills only ~1.5% more of the survivors near 2^65. Current state: 4 p at ~2^65 (~1 CPU-hour each) and 95 p at ~2^60. The next bit on the 95 is ~30× cheaper per kill than the next bit on the 4.
**Steps**
1. Priority = expected kills / expected cost = `(1/b) / 2^(b−1)·t_k` per p.
2. Pop the lowest-b p first (heap); push it back with b+1.
3. Replace the per-p deepening loop in `p_hunter_eliminate.py` (`process_p`, `levels`).

### C5. Cost/benefit stop rule
Stop TF when the marginal cost per kill exceeds the cost of a PRP test at that size (~hours on GPU, days on this CPU for p≈1.36e8). Compute with the measured ns/k rates. **Steps:** record seconds-per-bit in state; stop when `time_for_next_bit / (1/b) > time_per_prp_test`.

### C6. Complementary tests instead of ever-deeper TF
After a sensible TF depth, P-1 (Prime95 / `gwnum`) finds factors that TF can't reach cheaply, and a PRP test (Prime95 on CPU, GpuOwl on GPU) is the only thing that can actually *prove* a candidate. Export survivors as Prime95 `worktodo.txt` lines (`Pfactor=...`, `PRP=...`).

---

## Tier D: External bypass / hardware

1. **GPU TF (mfaktc/mfakto)**: Orders of magnitude faster than CPU TF (a mid-range GPU is ≈ 100–1000× this CPU). This machine has only an Intel UHD iGPU, so: (a) use `mfakto` on OpenCL if the iGPU is supported (expect small gain), or (b) run on a rented/other-machine NVIDIA GPU with `mfaktc` using a worktodo line `Factor=P,70,77`. Steps: export the surviving p list to worktodo format, run, import results.
2. **Custom wgpu/OpenCL kernel** for the iGPU: possible but 64×64→128 multiplication is emulated on Intel GPUs; expected gain is small versus A2 on 4 CPU cores. Low priority.
3. **Join PrimeNet**: results are credited and, importantly, your TF isn't duplicated.

---

## Tier E: Data and orchestration

1. **Persistent state in SQLite/JSONL**: `write_outputs` in `p_hunter_farm.py` rewrites the full JSON per batch (O(N) per batch, O(N²) overall). Use an append-only JSONL log or SQLite with one row per p.
2. **Rust/numpy segmented prime sieve** to replace the Python `next_primes` in `p_hunter_farm.py`.
3. **Store coverage as bit levels** (C1) so merging runs is `max()`.
4. **Add a `bench.py`** that runs fixed bands (53/58/64/71-bit) and prints ns/k so every change is measured the way CONTEXT.md demands.
5. **Extend `crosscheck.py`** for: Montgomery vs naive, q>2^64, known-factor regression (every factor in a pinned list must be recovered), and known Mersenne primes never being eliminated.

---

## Tier F: Considered and rejected

| Idea | Why it doesn't help |
|---|---|
| Cross-p shared sieving / q-major enumeration | The number of (p,k) pairs is the same; the sieve work per pair is unchanged. |
| Product-modulus batch modexp + GCD | Big-number FFT overhead is larger than constant-size mulmod. |
| Quartic/cubic power-residue pre-tests | Costs about as much as the modexp it'd avoid. |
| Smooth-k only | Success rate is ≈ 1/(2k) regardless of k's smoothness, so no hit-rate gain. |
| Sieving to √q | Removes only ~33% more modexps at heavy sieve cost. |
| Deeper sieving (>20M) now | 66% of survivors are prime; little left to remove (revisit only if A1/A2 don't land). |

---

## Prioritized roadmap

| # | Item | Impact | Effort |
|---|---|---|---|
| 1 | C2 — import PrimeNet TF/composite status, drop resolved p | **Huge** (may erase most of the work) | Low |
| 2 | A5 — `target-cpu=native`, LTO | Small–medium | Trivial |
| 3 | A1 — Montgomery u64 + A3 skip squarings | ~5–10× below 2^64 | Medium |
| 4 | A2 — 2-limb Montgomery | ~8–20× above 2^64 (removes the cliff) | Medium–High |
| 5 | C1 + C4 — bit-level ledger, breadth-first | Large efficiency | Medium |
| 6 | B1 + B2 — wheel + cache-sized bitset | 3–5× on sieve | Medium–High |
| 7 | A4 — interleaved modexp | 1.5–2× | Medium |
| 8 | B4/B5/B6/B7 — tuning, session cache, GIL, retire `p_hunter_rs` | 10–30% | Low–Medium |
| 9 | Tier D — GPU TF via mfaktc, P-1/PRP | 100×+ (on proper hardware) | Needs hardware |
| 10 | Tier E — SQLite, bench, crosscheck | Quality/robustness | Low |

**Rough combined outcome:** A1+A2+A5 plus B1/B2 should plausibly bring the deep regime (~200 ns/k) down to roughly 10–20 ns/k per core, i.e. 10–20×. Starting at the already-known TF depth (C2) multiplies that by avoiding the redundant bit levels entirely.

## Honest odds
GIMPS has already covered this range heavily (TF to ~2^77 for many exponents, and many are PRP-tested composite). The chance that a random unresolved exponent near 1.36e8 is a Mersenne prime is about `e^γ·ln(2^b)/p`-scale — on the order of 1 in tens of thousands per candidate even after TF, and your 100 survivors have mostly been (or will quickly be) resolved by others. The project is valuable as an engineering/learning exercise and for finding factors; treat a discovery as a lottery ticket and measure success by **ns/k**, **bit levels covered beyond the known frontier**, and **zero missed known factors**.
