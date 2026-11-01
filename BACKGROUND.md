# Background: what is `k`, and how does this map to GIMPS terms?

You don't need a math degree to use this project, but a few ideas make every
menu, option and report make sense. This page explains them plainly, then maps
the words used here to the words used by [GIMPS](https://www.mersenne.org/)
(the Great Internet Mersenne Prime Search) and its tools.

## 1. The goal

A **Mersenne number** is `M_p = 2^p - 1`. It can only be prime if the exponent
`p` is itself prime, but most prime `p` still give composite (non-prime)
numbers. For example `M_11 = 2047 = 23 × 89`. Only a few dozen Mersenne primes
are known. This project tries to *throw away* exponents that are provably
composite, as cheaply as possible, so that expensive primality tests are only
run on the few that are left.

The cheapest way to prove `M_p` composite is to **find one factor**. That is
called **trial factoring (TF)**.

## 2. What `k` is

Any prime factor `q` of `2^p - 1` (with `p` prime) must have the form:

```text
q = 2 · k · p + 1        for some whole number k = 1, 2, 3, ...
```

So instead of testing every number as a possible factor, you only test
`q = 2p+1, 4p+1, 6p+1, ...`. **`k` is simply the counter that walks along
those candidates.** `k = 1` is the first candidate, `k = 2` the second, and so
on.

Why this form? If `q` divides `2^p - 1`, then `2^p ≡ 1 (mod q)`. The
"order" of 2 modulo `q` must divide `p`, and since `p` is prime and the order
isn't 1, it equals `p`. By Fermat's little theorem the order also divides
`q - 1`, so `p` divides `q - 1`. Because `q - 1` is even, `2p` divides it, which
gives `q = 2kp + 1`.

Worked examples:

| Mersenne number | Factor `q` | `q = 2kp+1` | `k` |
| --- | --- | --- | --- |
| `M_11 = 2047` | 23 | `2·1·11 + 1` | 1 |
| `M_11 = 2047` | 89 | `2·4·11 + 1` | 4 |
| `M_23` | 47 | `2·1·23 + 1` | 1 |
| `M_29` | 233 | `2·4·29 + 1` | 4 |

### The two extra rules that shrink the search

1. **`q` must be ≡ 1 or 7 (mod 8).** (2 must be a quadratic residue modulo
   `q`.) This throws out about half of all `k` immediately.
2. **`q` must be prime**, so any `k` where `q` is divisible by a small prime
   (3, 5, 7, 11, 13, ...) is skipped. This is the **sieve**. In this project
   the sieve walks primes up to about 20 million, and only the `k` that survive
   get the real test.

### The real test (the "modexp")

For each surviving `k`, compute `2^p mod q`. If the answer is `1`, then `q`
divides `M_p`, so you found a factor and `M_p` is composite. This is a
modular exponentiation done by repeated squaring (about 27 squarings for
`p ≈ 1.36e8`). Making this step fast is the whole game; see
[OPTIMIZATION_IDEAS.md](OPTIMIZATION_IDEAS.md).

## 3. How big does `k` get? ("depth" and bit levels)

"Depth" here means "how far along `k` have we searched?" There are two ways to
say it:

- **As a `k` limit**: "tested every `k ≤ 1,000,000,000`" (this project's
  `--depth-k` and `max_k` options).
- **As a bit level** (the GIMPS way): "no factor below `2^b`."

They are the same thing in different units, because `q ≈ 2kp`:

```text
bit level b  ≈  log2(2 · k · p)          k_max ≈ 2^b / (2p)
```

For `p = 136,279,843`:

| `k` limit | Largest `q` tested | Bit level |
| --- | --- | --- |
| 1 × 10⁹ | about 2.7 × 10¹⁷ | about 2^58 |
| 1 × 10¹¹ | about 2.7 × 10¹⁹ | about 2^64.5 |
| 1 × 10¹³ | about 2.7 × 10²¹ | about 2^71 |
| 2^77 level | about 1.5 × 10²³ | `k` ≈ 5.5 × 10¹⁴ |

**Each extra bit doubles the amount of work** but finds only a few percent more
factors, which is why deep trial factoring has diminishing returns. This is
also why this project's speed drops when `q` passes `2^64`: the numbers no
longer fit in a single 64-bit register.

## 4. Glossary: this project ↔ GIMPS

| Term here | GIMPS / community term | Meaning |
| --- | --- | --- |
| `p` | exponent | The `p` in `2^p - 1`. |
| `M_p` | Mersenne number | `2^p - 1`. |
| `k`, `max_k`, `--depth-k` | (implicit in "bit level") | The candidate counter in `q = 2kp + 1`. |
| factor probe, `tf_window`, "gate" | trial factoring (TF) | Testing candidate factors `q`. |
| depth in `k` | TF bit depth ("TF to 2^70") | How far TF has searched. |
| survivor / passed | "no factor found", unfactored | Not eliminated *yet*. Not a proof of primality. |
| dead / eliminated | factored | A factor was found; `M_p` is composite. |
| sieve | TF sieve (mfaktc "sieve primes") | Removes `k` whose `q` has a small prime factor. |
| deepening | extending TF to a higher bit level | Searching more `k` on survivors. |
| farm | (none) | Bulk scan of many exponents to a fixed depth. |
| eliminator | (none) | Spends extra TF time per surviving exponent. |
| close-out primality check | PRP / LL test | Final primality test. Here only a quick check; use Prime95/GpuOwl for real. |
| likeness score | (none) | A heuristic ranking from the binary shape. Not a proof and not a filter. |
| known exponent | known Mersenne prime exponent | One of the ~52 exponents proven prime. Used as a regression test. |

Other GIMPS words you will run into elsewhere:

- **PrimeNet**: the server that hands out and records GIMPS work.
- **Prime95 / mprime**: the standard CPU software (TF, P-1, PRP/LL).
- **mfaktc / mfakto**: GPU trial-factoring programs (NVIDIA / AMD-OpenCL).
- **GpuOwl**: GPU PRP testing.
- **P-1 factoring**: another factor-finding method that finds factors TF
  can't reach cheaply. Complementary to TF.
- **LL (Lucas-Lehmer) test**: the classic *proof* of primality for `M_p`.
- **PRP (probable prime) test**: today's faster standard first-pass test; a
  Mersenne PRP result is then confirmed.
- **GHz-days**: GIMPS's unit for measuring work done.
- **TF target / "TF to 2^77"**: how deep TF is considered worthwhile for an
  exponent. For `p` around 1.36 × 10⁸, PrimeNet's target is about `2^77`.

## 5. Reading a result

- *"p = X survived to k ≤ 1e9"* → nobody found a factor with `q` up to about
  `2·1e9·X`. Many such exponents are still composite; the factor is just
  larger.
- *"p = X dead, factor q = ..."* → `M_X` is composite. Done.
- *"p = X passed (eliminator)"* → survived this project's chosen time/depth
  budget. Still not a proof.

## 6. Where this sits in the bigger GIMPS pipeline

```text
Choose exponent p  →  trial factoring (TF)  →  P-1 factoring  →  PRP / LL test
   (cheap per p)       (this project)           (Prime95)         (days of CPU/GPU)
```

GIMPS has already pushed TF very deep on most exponents, mostly with GPUs. The
practical value of a CPU project like this is to **find new factors, test
ideas for going faster, and check which exponents still have shallow TF** (see
section C2 of [OPTIMIZATION_IDEAS.md](OPTIMIZATION_IDEAS.md)). It is worth
checking an exponent's status on mersenne.org before investing time in it.
