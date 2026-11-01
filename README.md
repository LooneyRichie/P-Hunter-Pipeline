# p_hunter_pipeline

`p_hunter_pipeline` is a collection of Python command-line tools and optional
Rust/PyO3 extensions for searching prime exponents `p` whose Mersenne numbers

```text
M_p = 2^p - 1
```

may be prime. It screens exponents, searches for small factors, deepens
trial-factoring coverage, checkpoints long-running work, and reports survivors
for later primality testing.

## The story

I've chased Mersenne primes for more than ten years. This project started as
a small Python script I carried over from my phone with the bare minimum of
tooling, and it grew into a sieve, a farm, an eliminator and two Rust
trial-factoring engines. It is an unorthodox, hands-on hobby project built
around one belief: at this scale the normal approach stops working, so you have
to think differently about where the time goes. I'm publishing it in case
someone finds it useful or can make it better. Pull requests, ideas and
benchmarks are welcome.

## The challenge

Finding a new Mersenne prime sounds impossible, and for a lone person with a
laptop-class CPU the odds are long. But every big leap in this field started
with someone who saw a shortcut everyone else missed: a smarter sieve, a better
multiplication, a cleverer way to skip work. This project was built on a
4-core mini PC, not a supercomputer, and it already found ways to throw away
huge amounts of wasted effort. You don't need a supercomputer to change the
world. You need a good idea.

**So here is the invitation.** If you love problems that look impossible, take
a swing:

- Find a faster way through the slow parts. The biggest known bottleneck is
  the jump in cost above 2^64 (see [Performance snapshot](#performance-snapshot)).
- Find work that can be skipped entirely: a smarter starting depth, a better
  list of exponents, or a pattern nobody has used yet.
- Beat the numbers in [`OPTIMIZATION_IDEAS.md`](OPTIMIZATION_IDEAS.md), or
  prove one of its ideas wrong and tell us why.
- Open an issue or pull request with your benchmark. A small, measured
  improvement counts. The seed of the next big idea might be yours.

Run it, break it, improve it. If you find something, credit is requested in
the [LICENSE](LICENSE), and the discovery is yours to celebrate.

> **A survivor is not a proof.** “No factor found” means only that the selected
> trial-factor search did not find a factor in the searched range. Use a
> Lucas-Lehmer test (for example, through GIMPS/Prime95) to establish whether
> `M_p` is prime.

> **New here?** [INSTRUCTIONS.md](INSTRUCTIONS.md) tells you exactly which
> script to run for what, and [BACKGROUND.md](BACKGROUND.md) explains what `k`
> is, how depth and bit levels work, and maps this project's terms to GIMPS
> terms. Install Python dependencies with
> `python -m pip install -r requirements.txt`.

## Contents

- [The story](#the-story)
- [The challenge](#the-challenge)
- [How the search works](#how-the-search-works)
- [Requirements and setup](#requirements-and-setup)
- [Choose a workflow](#choose-a-workflow)
- [Commands and examples](#commands-and-examples)
- [Files and checkpoints](#files-and-checkpoints)
- [Usage tips](#usage-tips)
- [Limitations and correctness notes](#limitations-and-correctness-notes)
- [Performance snapshot](#performance-snapshot)
- [Future improvements](#future-improvements)
- [Contributing](#contributing)
- [License and credit](#license-and-credit)

## How the search works

For an odd prime exponent `p`, any prime factor `q` of `M_p` must have the form
`q = 2kp + 1` for some positive integer `k`. The trial-factoring routines
search for such `q` values and test whether `2^p ≡ 1 (mod q)`.

The main candidate sieve in [`revised_binV3.py`](revised_binV3.py) applies
necessary checks, tests that `p` is prime, probes for factors through
`SCAN_MAX_K` (currently `100,000,000`), and applies a last-digit check.
Binary-shape likeness is calculated separately as a **ranking heuristic**; it
is not a proof and does not filter candidates in this pipeline.

The larger-range workflow is split into stages:

1. [`p_hunter_farm.py`](p_hunter_farm.py) scans prime exponents and trial-factors
   each to a chosen depth.
2. [`p_hunter_range.py`](p_hunter_range.py) scans a specified interval and
   merges its survivors into the farm survivor list.
3. [`p_hunter_eliminate.py`](p_hunter_eliminate.py) spends more factoring time
   on survivors, recording candidates eliminated by a factor and those that
   passed the current factoring budget.

The interactive sieve, iterative deepener, and browser interface provide
alternative workflows for smaller or exploratory scans.

## Requirements and setup

- Python 3.8 or newer is declared by the `p_hunter_tf` Python build
  configuration.
- The recommended high-throughput factor-probe backend is
  [`p_hunter_tf/`](p_hunter_tf/), a Rust extension built with Cargo and
  [maturin](https://www.maturin.rs/).
- `p_hunter_rs/` is an older optional Rust backend. Its factor values are
  limited to 64 bits, and some callers explicitly cap its usable depth.
- The Python-only paths are available in some tools, but can be much slower.
- `scikit-learn` is only needed to fit optional scorer calibration weights.
  Close-out primality checks try `gmpy2`, `numthy`, or `sympy` when available.

Run commands from the project root so imports and the default state/output
paths resolve consistently. This workspace already has a `.venv`; use it rather
than creating a second environment:

```bash
source .venv/bin/activate
python -m pip install "maturin>=1,<2"
cd p_hunter_tf
maturin develop --release
cd ..
```

The maturin step builds and installs the `p_hunter_tf` extension into the
active Python environment. It requires a working Rust/Cargo toolchain.
`p_hunter_tf` is recommended for the farm, range, and eliminator workflows.

To build the older backend instead, activate the same environment and run:

```bash
cd p_hunter_rs
maturin develop --release
cd ..
```

The backend is optional for the interactive candidate sieve and
`p_hunter_deepen.py`, which have Python fallbacks. The farm requires a Rust
backend; the eliminator requires `p_hunter_tf`.

## Choose a workflow

| Goal | Entry point | Notes |
| --- | --- | --- |
| Explore sieve features, run a benchmark, or inspect one exponent | [`revised_binV3.py`](revised_binV3.py) | Interactive menu; includes range, endless, parallel, and deep-verify options. |
| Iteratively factor a small range | [`p_hunter_deepen.py`](p_hunter_deepen.py) | Command-line options; can use `p_hunter_rs`, otherwise falls back to Python. |
| Scan above the largest known exponent for a fixed number of survivors | [`p_hunter_farm.py`](p_hunter_farm.py) | Multiprocess farm with a resumable checkpoint. |
| Scan a specific high-exponent interval | [`p_hunter_range.py`](p_hunter_range.py) | Saves a per-range checkpoint and merges results into the farm survivor JSON. |
| Spend a deeper, bounded factoring budget on farm survivors | [`p_hunter_eliminate.py`](p_hunter_eliminate.py) | Interactive menu; requires `p_hunter_tf`. |
| Use the browser interface to configure a deepening scan | [`p_hunter_web.py`](p_hunter_web.py) | Local HTTP server; results can be downloaded as JSON. |

`base/` contains earlier snapshots of the deepening script, eliminator, and
candidate sieve. The top-level files are the current entry points; use the
`base/` copies only when comparing historical behavior.

## Commands and examples

All examples below assume the project root is the current directory and the
existing virtual environment is active. They show the current command-line
interfaces; check `--help` on scripts that provide it before starting a long
run.

### Interactive candidate sieve and benchmark

```bash
python revised_binV3.py
```

Choose from the menu to run a range search, endless search, resume a saved
search, benchmark against the known-exponent list, analyze one `p`, run a
checkpointed deep verification, or search a range in parallel. The benchmark
reports known exponents captured, unmatched survivors, and elapsed time.

### Iterative deepening over a short range

```bash
python p_hunter_deepen.py \
  --start 10000 --end 10050 \
  --k 1000 --mult 10 --rounds 3 \
  --workers 4 --save results/deepen.json
```

Each round tests the next `k` band on remaining candidates. `--interactive`
prompts for the next maximum `k` instead of using `--mult`. The script also
asks whether to run a close-out primality test after reporting; answer `n` to
skip it.

### Farm search with a Rust backend

```bash
python p_hunter_farm.py \
  --target 100 --depth-k 1000000000 --workers 4
```

The defaults start above `136,279,841`, seek 100 new survivors, and factor
through `k <= 1,000,000,000`. To continue an existing farm checkpoint:

```bash
python p_hunter_farm.py --resume --target 100 --workers 4
```

The farm stores progress in `farm_state.json` by default and writes
`farm_survivors.json` and `farm_survivors.txt`. `--state PATH` selects a
different checkpoint path, but the output filenames remain fixed in the
current working directory.

### Scan and merge a particular interval

```bash
python p_hunter_range.py 150000000 150050000 100000000
```

The positional arguments are start exponent, end exponent, and optional
maximum `k`. This scan creates a `range_<start>_<end>.json` checkpoint and
merges results into `farm_survivors.json`; existing survivors are retained.
The range tool also supports an interactive mode when called without
positional arguments.

### Deepen the farm survivors

```bash
python p_hunter_eliminate.py
```

Use menu option **1** to start or continue factoring, **2** to view the list,
**3** to adjust the per-exponent time limit, starting depth, ceiling, and
worker count, and **4** to re-queue entries marked passed. The eliminator
imports new entries from `farm_survivors.json` and updates
`eliminate_state.json` plus the deepened JSON and text reports.

### Browser interface

```bash
python p_hunter_web.py
```

By default it serves the UI on `http://127.0.0.1:8000` and attempts to open a
browser. To choose a port and avoid opening a browser automatically:

```bash
python p_hunter_web.py --port 9000 --no-browser
```

The page configures the exponent range, starting `k`, multiplier, and worker
count. After a run, it can show the survivors, start a close-out primality
check, and download the run status/results as JSON.

## Files and checkpoints

These scripts use paths relative to the current working directory. Keep a
backup of any state/results you want to preserve before starting a run that
may update them.

| File | Purpose |
| --- | --- |
| `farm_state.json` | Farm checkpoint: next exponent, accumulated survivors, and scan counters. |
| `farm_survivors.json` | Shared farm/range survivor list consumed by the eliminator. |
| `farm_survivors.txt` | Human-readable farm survivor report. |
| `range_<start>_<end>.json` | Checkpoint for a specific interval scanned by `p_hunter_range.py`. |
| `eliminate_state.json` | Per-exponent factoring coverage, status, elapsed work, and settings. |
| `farm_survivors_deepened.json` | Deepened records grouped as `pending`, `passed`, and `dead`. |
| `farm_survivors_deepened.txt` | Human-readable deepened report. |
| `p_hunter_state.json` | Checkpoint for the endless search in `revised_binV3.py`. |
| `deep_verify_state.json` | Checkpoint for the interactive sieve's deep-verify operation. |
| `scorer_calibration.json` | Optional weights generated by the scorer calibration menu option. |

The farm's `--resume` option resumes its selected state file; it does not start
a fresh scan. The eliminator's state is cumulative, and the range scanner
merges results rather than replacing the survivor list. Do not remove or
overwrite checkpoint files unless you intend to discard that progress.

## Usage tips

- **Start with a short range.** Use a narrow exponent range and low `k` depth
  to verify backend installation and understand throughput before scaling up.
- **Set worker counts deliberately.** Farm and interval scans default to the
  available CPU count. Reducing `--workers` can leave capacity for other
  processes and reduce memory pressure.
- **Factor in stages.** A shallow farm/range scan is a prefilter; the
  eliminator is designed to spend additional time on the smaller survivor
  list. Increasing depth costs more work per exponent.
- **Resume matching work.** Keep the relevant checkpoint and rerun with the
  same range/configuration where applicable. The web UI is in-memory and does
  not provide the farm's persistent checkpoint behavior.
- **Use the faster backend for deep factoring.** Build `p_hunter_tf` before
  using the farm, range scanner, or eliminator for substantial work. If the
  farm falls back to `p_hunter_rs`, it warns that the usable `k` range is
  capped.
- **Treat likeness scores as prioritization only.** A higher score does not
  establish primality, and a low score is not grounds to reject a survivor.
- **Verify promising results independently.** Use a Lucas-Lehmer test or
  trusted Mersenne-primality software for a final result.
- **Keep runs local.** The web server defaults to loopback. Avoid binding it
  to a public or untrusted network; the interface has no authentication.

## Limitations and correctness notes

- A “passed” or “survived” record means no factor was found through the
  recorded `k` depth. It does not mean `M_p` is prime.
- The core candidate sieve's scan depth is set by `SCAN_MAX_K` in
  `revised_binV3.py`; farm and interval scan depths are controlled separately.
- The `p_hunter_tf` extension uses a segmented `k` window and supports factor
  candidates below `2^100`; a requested window whose largest `q` reaches
  `2^100` is rejected.
- The older `p_hunter_rs` extension uses 64-bit factor values and should not
  be treated as interchangeable with the deeper `p_hunter_tf` backend.
- `p_hunter_deepen.py` can run with its Python fallback, but that path may be
  impractical for large ranges or high factoring depths.
- The optional close-out checks select an available primality-test library.
  For mathematical confirmation of large Mersenne numbers, use a
  Lucas-Lehmer implementation rather than interpreting “passed” as proof.
- There is no top-level dependency manifest or unified automated test suite.
  The `p_hunter_tf/crosscheck.py` script cross-checks the Rust factor window
  against brute force on selected ranges and known cases.

## Performance snapshot

Single-core `p_hunter_tf` `tf_window` timings at `p = 136279843` on an Intel
N200 (4 cores, AVX2, no AVX-512), 5e7-wide `k` bands:

| `k` band | Candidate size `q` | ns per `k` |
| --- | --- | --- |
| up to 5e7 | about 2^53 | 11.3 (a factor was found early) |
| about 1e9 | about 2^58 | 41.9 |
| about 5e10 | below 2^64 | 43.2 |
| about 1e11 | above 2^64 | 201.9 |
| about 1e13 | about 2^71 | 210.8 |

Crossing 2^64 makes each `k` about 4.7 times slower because the wide
modular multiplication falls back to software u128 division. Fixing that
(Montgomery arithmetic) is the biggest single win available.

**Check what is already known.** PrimeNet reports exponents in this range, for
example `136282409`, as already proven composite with trial factoring done to
2^77. Before spending CPU time on an exponent, check its status at
[mersenne.org](https://www.mersenne.org/) so you only work on bit levels that
nobody has covered. Treat any new Mersenne prime as a lottery ticket; the
value here is the engineering and the factors found along the way.

The full list of speedup ideas, with steps, validation and a prioritized
roadmap, is in [`OPTIMIZATION_IDEAS.md`](OPTIMIZATION_IDEAS.md).

## Future improvements

0. **Implement the roadmap in [`OPTIMIZATION_IDEAS.md`](OPTIMIZATION_IDEAS.md)**:
   import known trial-factoring depth from PrimeNet, Montgomery arithmetic
   (64-bit and multi-limb), a k-class wheel with a cache-sized bitset sieve,
   a per-exponent bit-level ledger with breadth-first scheduling, and GPU
   trial factoring or P-1/PRP hand-off.
1. **Add automated regression tests** for gate invariants, known Mersenne
   exponents, checkpoint/resume behavior, and factor-window boundaries. Make
   zero missed known exponents a required correctness condition.
2. **Standardize benchmarks** so candidate recall, survivor counts, factoring
   depth, elapsed time, and throughput can be compared consistently before
   and after sieve changes.
3. **Share the core candidate and factor-probe logic** among the top-level
   scripts, and explicitly label or retire historical copies under `base/` to
   reduce behavioral drift.
4. **Cache larger prime sieves efficiently**, including persistent reuse where
   worthwhile, to reduce repeated setup cost in long and multi-process runs.
5. **Harden checkpoint management** with explicit output directories,
   configuration/schema validation, and clearer safeguards against resuming
   a checkpoint with incompatible settings.
6. **Improve final verification integration** by clearly separating
   probabilistic primality checks from proof-producing Lucas-Lehmer testing
   and by recording verification provenance with each result.
7. **Make backend selection and limits uniform** across the farm, range,
   deepening, and interactive tools, with a small backend conformance suite.
8. **Expand operational reporting** with per-stage timings, a concise run
   summary, and export formats that preserve the searched factor bounds and
   backend version.
