//! Minimal PyO3 bridge: sieved B2 probe with min_k support.

use std::sync::{Mutex, OnceLock};
use pyo3::prelude::*;

// ---------- helpers ----------

fn isqrt_u128(n: u128) -> u128 {
    if n < 2 { return n; }
    let mut x = (n as f64).sqrt() as u128;
    while x.checked_mul(x).map_or(true, |xx| xx > n) { x -= 1; }
    while (x + 1).checked_mul(x + 1).map_or(false, |xx| xx <= n) { x += 1; }
    x
}

fn mod_pow_u64(mut base: u64, mut exp: u64, m: u64) -> u64 {
    if m == 1 { return 0; }
    let mut r: u64 = 1;
    base %= m;
    while exp > 0 {
        if exp & 1 == 1 { r = ((r as u128 * base as u128) % m as u128) as u64; }
        exp >>= 1;
        base = ((base as u128 * base as u128) % m as u128) as u64;
    }
    r
}

fn mod_pow_u128(mut base: u128, mut exp: u128, m: u128) -> u128 {
    if m == 1 { return 0; }
    let mut r: u128 = 1;
    base %= m;
    while exp > 0 {
        if exp & 1 == 1 { r = (r * base) % m; }
        exp >>= 1;
        base = (base * base) % m;
    }
    r
}

fn sieve_primes(limit: u64) -> Vec<u32> {
    if limit < 2 { return Vec::new(); }
    let mut s = vec![true; (limit + 1) as usize];
    s[0] = false; s[1] = false;
    let mut i = 2u64;
    while i * i <= limit {
        if s[i as usize] {
            let mut j = i * i;
            while j <= limit { s[j as usize] = false; j += i; }
        }
        i += 1;
    }
    s.iter().enumerate()
        .filter_map(|(i, &b)| (b && i >= 2).then_some(i as u32))
        .collect()
}

static PRIME_CACHE: OnceLock<Mutex<Vec<u32>>> = OnceLock::new();

fn cached_primes(limit: u64) -> Vec<u32> {
    let cell = PRIME_CACHE.get_or_init(|| Mutex::new(Vec::new()));
    let mut cache = cell.lock().unwrap();
    let have = cache.last().copied().unwrap_or(0) as u64;
    if have < limit {
        *cache = sieve_primes(limit);
    }
    let n = cache.partition_point(|&x| (x as u64) <= limit);
    cache[..n].to_vec()
}

// ---------- the exported function ----------

#[pyfunction]
#[pyo3(signature = (p, min_k, max_k))]
fn gate_b2_factor_probe(p: u64, min_k: u64, max_k: u64) -> Option<(u64, u64)> {
    if p < 8 || p & 1 == 0 { return None; }
    if min_k > max_k { return None; }

    let k_max = max_k;
    let q_max: u128 = 2 * (k_max as u128) * (p as u128) + 1;

    let mut alive = vec![1u8; (k_max as usize) + 1];
    alive[0] = 0;

    let s_limit = isqrt_u128(q_max);
    const SIEVE_CAP: u64 = 50_000_000;
    let sieve_limit = s_limit.min(SIEVE_CAP as u128) as u64;

    let primes = cached_primes(sieve_limit);
    let two_p = 2u128 * (p as u128);

    for &s_u32 in &primes {
        let s = s_u32 as u128;

        // q = 2kp+1 is always odd, so s=2 has no root in k.  Without
        // this skip, the inverse call below degenerates and marks all
        // odd k as composite.
        if s == 2 { continue; }

        // q ≡ 1 (mod p) always: no root in k.
        if s == p as u128 { continue; }

        let inv = mod_pow_u128(two_p % s, s - 2, s);
        let k0  = s - inv;                       // in [1, s-1]

        // If q(k0) == s, then k0 marks a *prime* q, not a composite.
        // Skip one period so q = s itself survives the sieve.
        let q0 = 2 * k0 * (p as u128) + 1;
        let first = if q0 == s { k0 + s } else { k0 };

        if first > k_max as u128 { continue; }
        for slot in alive[(first as usize)..].iter_mut().step_by(s as usize) {
            *slot = 0;
        }
    }

    let start = min_k.max(1);
    for k in start..=k_max {
        if alive[k as usize] == 0 { continue; }
        let q128 = 2u128 * (k as u128) * (p as u128) + 1;

        if (q128 & 7) != 1 && (q128 & 7) != 7 { continue; }
        if q128 > u64::MAX as u128 { continue; }

        let q = q128 as u64;
        if mod_pow_u64(2, p, q) == 1 {
            return Some((q, k));
        }
    }
    None
}

#[pymodule]
fn p_hunter_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(gate_b2_factor_probe, m)?)?;
    Ok(())
}