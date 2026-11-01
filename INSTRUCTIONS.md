# Which script do I run?

A practical guide to using p_hunter_pipeline. For background and the math see
[README.md](README.md); for speedup ideas see
[OPTIMIZATION_IDEAS.md](OPTIMIZATION_IDEAS.md). New to the terms (`k`, bit
level, TF, PRP)? Read [BACKGROUND.md](BACKGROUND.md) first.

## 0. One-time setup

Run everything from the project root, inside a virtual environment.

```bash
python -m venv .venv            # skip if .venv already exists
source .venv/bin/activate       # Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt

# Build the fast Rust trial-factoring engine (needs Rust: https://rustup.rs)
cd p_hunter_tf
maturin develop --release
cd ..

# Optional older 64-bit engine (only the deepen script and fallbacks use it)
cd p_hunter_rs
maturin develop --release
cd ..
```

Check that it worked:

```bash
python -c "import p_hunter_tf; print('p_hunter_tf OK')"
```

Rust is optional for some tools: `revised_binV3.py` and `p_hunter_deepen.py`
fall back to pure Python (slower). The farm needs a Rust backend and the
eliminator needs `p_hunter_tf`.

Verify the Rust kernel any time with:

```bash
python p_hunter_tf/crosscheck.py
```

## 1. Pick your tool

| I want to... | Run | Rust needed? | Saves progress? |
| --- | --- | --- | --- |
| Explore with a menu, benchmark, check one exponent | `python revised_binV3.py` | Optional (Python fallback) | Yes (menu options) |
| Get a nicely formatted HTML view | `python p_hunter_web.py` | Optional | No (in memory, JSON download) |
| Grind a list of exponents with growing `k` until few remain | `python p_hunter_deepen.py` | Optional | Only with `--save` |
| Scan a huge range for survivors, save them to a file | `python p_hunter_farm.py` | Yes | Yes (`farm_state.json`) |
| Scan one specific exponent interval and merge results | `python p_hunter_range.py` | Yes | Yes (per-range file) |
| Spend serious time deepening the farm's survivors | `python p_hunter_eliminate.py` | Yes (`p_hunter_tf`) | Yes (`eliminate_state.json`) |

Typical big-run workflow: **farm → (range) → eliminate**. Beginners should
start with `revised_binV3.py` or the web UI.

---

## 2. `revised_binV3.py` — the menu system

The core engine and the friendliest entry point. It uses the Rust backend when
installed and quietly falls back to Python when not.

```bash
python revised_binV3.py
```

What the menu offers: range search, endless search, resume a saved search,
a benchmark against the known Mersenne exponents, analysis of a single `p`,
checkpointed deep verification, and a parallel range search.

Tips:
- Run the **benchmark** first. It reports how many known exponents are
  captured; a healthy build captures all of them.
- Start with a small range to see how fast your machine is.
- Scores ("likeness") only rank survivors. They never prove anything.
- Other scripts import this file, so keep it in the project root.

## 3. `p_hunter_web.py` — HTML dashboard

```bash
python p_hunter_web.py                           # http://127.0.0.1:8000
python p_hunter_web.py --port 9000 --no-browser  # pick a port, no auto-open
```

Set the exponent range, starting `k`, multiplier and worker count in the page,
run, then view the formatted survivors and download the results as JSON.

Know the drawbacks:
- The close-out primality check needs **`gmpy2`** for huge exponents. Without
  it the page falls back to `numthy`/`sympy`, which are far too slow for large
  `p`. Install it with `pip install gmpy2`.
- State is **in memory only**: download your JSON before closing the server.
- It has no authentication. Keep it on `127.0.0.1`; do not expose it to a
  network.

## 4. `p_hunter_deepen.py` — iterative deepening on a list of exponents

Runs a set of `k` values over a range of exponents, then multiplies `k` and
repeats on whatever is left, until nothing more is eliminated or few candidates
remain.

```bash
python p_hunter_deepen.py --start 10000 --end 10050 \
    --k 1000 --mult 10 --rounds 3 --workers 4 --save results/deepen.json
```

| Option | Meaning |
| --- | --- |
| `--start`, `--end` | Exponent range (default 10000–10050) |
| `--k` | Starting maximum `k` (default 1000) |
| `--mult` | Multiplier applied to `k` each round (default 10) |
| `--rounds` | Stop after this many rounds (default: until it stalls) |
| `--interactive` | Ask for the next maximum `k` each round instead of `--mult` |
| `--workers` | Parallel processes (default: auto) |
| `--save` | Write results to this JSON file (create the folder first) |

Afterwards it offers a close-out primality test; answer `n` to skip it.
Best for small and medium ranges and for learning how survivor counts fall as
`k` grows.

## 5. `p_hunter_farm.py` — big-range scan, saved for later

Scans prime exponents upward from a start point and trial-factors each to a
fixed depth, saving survivors for the eliminator.

```bash
python p_hunter_farm.py --target 100 --depth-k 1000000000 --workers 4
python p_hunter_farm.py --resume --target 100 --workers 4   # continue later
```

| Option | Meaning |
| --- | --- |
| `--start` | First exponent (default: just above 136,279,841) |
| `--target` | New survivors to collect this run (default 100) |
| `--depth-k` | Trial-factor every `k` up to this (default 1e9) |
| `--workers` | Processes (default: all CPUs; use physical cores) |
| `--chunk` | Primes per batch (default 4 × workers) |
| `--resume` | Continue from the checkpoint instead of starting fresh |
| `--state` | Checkpoint path (default `farm_state.json`) |

Outputs: `farm_state.json` (checkpoint), `farm_survivors.json` and
`farm_survivors.txt` (the survivor list, created in the current directory).

Always pass `--resume` when continuing. Without it you restart from the
beginning.

## 6. `p_hunter_range.py` — scan a specific interval

```bash
python p_hunter_range.py 150000000 150050000 100000000
#                        start     end       max k (optional)
```

Writes a `range_<start>_<end>.json` checkpoint and **merges** its survivors
into `farm_survivors.json` without deleting existing ones. Run it with no
arguments for an interactive prompt.

## 7. `p_hunter_eliminate.py` — spend a deep budget on survivors

Reads `farm_survivors.json`, imports new entries, and gives each exponent a
bounded amount of deep trial-factoring, recording those killed by a factor and
those that passed.

```bash
python p_hunter_eliminate.py
```

Menu:
1. **Start / continue** factoring
2. **View** the list (pending / passed / dead)
3. **Settings**: per-exponent time limit, starting depth, ceiling, workers
4. **Re-queue** entries marked passed (to push them deeper)

Outputs: `eliminate_state.json`, `farm_survivors_deepened.json` and `.txt`.
State is cumulative, so you can stop (Ctrl+C) and resume safely.

Cost warning: each extra bit of factoring depth roughly doubles the work for a
small number of extra kills. See [OPTIMIZATION_IDEAS.md](OPTIMIZATION_IDEAS.md)
section C4/C5.

---

## 8. Recipes

**Learn the tool (5 minutes):**
`python revised_binV3.py` → benchmark → try a small range.

**Quick visual exploration:**
`python p_hunter_web.py` → pick a range → run → download JSON.

**Long unattended hunt:**
```bash
python p_hunter_farm.py --target 100 --depth-k 1000000000 --workers 4
python p_hunter_eliminate.py        # menu option 1, leave it running
```

**Add another region to the same survivor list:**
```bash
python p_hunter_range.py 150000000 150050000 100000000
python p_hunter_eliminate.py
```

**After you have survivors:** a survivor is *not* a proof. Check it on
[mersenne.org](https://www.mersenne.org/) first (many exponents are already
factored or tested), then run a Lucas-Lehmer / PRP test with Prime95 or GpuOwl.

## 9. Good habits

- Back up `*.json` state files before experimenting; they hold your progress.
- Set `--workers` to your physical core count.
- Run from the project root so relative paths and imports work.
- Change one thing at a time and benchmark before and after.
- Never treat "passed" as "prime".

## 10. Troubleshooting

| Problem | Fix |
| --- | --- |
| `ModuleNotFoundError: p_hunter_tf` | Activate the venv and run `maturin develop --release` in `p_hunter_tf/` |
| Farm says no Rust backend / depth capped | Build `p_hunter_tf` (the older `p_hunter_rs` is 64-bit only) |
| `maturin` or `cargo` not found | Install Rust from rustup.rs, then `pip install maturin` |
| `gmpy2` fails to install | Install the GMP/MPFR/MPC development packages (e.g. `apt install libgmp-dev libmpfr-dev libmpc-dev`) or rely on the slower fallbacks |
| Web page loses results | State is in memory; download the JSON before stopping the server |
| Run restarted from scratch | You forgot `--resume` (farm) or changed the state path |
| `p_hunter_rs` build error about edition 2024 | Update Rust: `rustup update` |
