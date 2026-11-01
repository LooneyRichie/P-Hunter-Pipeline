//! Segmented trial-factoring window for Mersenne numbers: does q = 2kp+1 divide 2^p - 1?
//! Memory is bounded by BLOCK regardless of k range; q up to 2^MAX_Q_BITS is handled exactly.

use std::sync::{Arc, Mutex, OnceLock};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

const BLOCK: u64 = 1 << 23;
const SIEVE_CAP: u64 = 20_000_000;
const MAX_Q_BITS: u32 = 100;

fn isqrt_u128(n: u128) -> u128 {
    if n < 2 { return n; }
    let mut x = (n as f64).sqrt() as u128;
    while x.checked_mul(x).map_or(true, |xx| xx > n) { x -= 1; }
    while (x + 1).checked_mul(x + 1).map_or(false, |xx| xx <= n) { x += 1; }
    x
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

static PRIME_CACHE: OnceLock<Mutex<(u64, Arc<Vec<u32>>)>> = OnceLock::new();

fn cached_primes(limit: u64) -> (Arc<Vec<u32>>, usize) {
    let cell = PRIME_CACHE.get_or_init(|| Mutex::new((0, Arc::new(Vec::new()))));
    let mut g = cell.lock().unwrap();
    if g.0 < limit {
        *g = (limit, Arc::new(sieve_primes(limit)));
    }
    let v = g.1.clone();
    let n = v.partition_point(|&x| (x as u64) <= limit);
    (v, n)
}

// inverse of a mod m (m prime, a % m != 0), result in [1, m-1]
fn inv_mod(a: u64, m: u64) -> u64 {
    let (mut r0, mut r1) = (m as i64, (a % m) as i64);
    let (mut t0, mut t1) = (0i64, 1i64);
    while r1 != 0 {
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (t0, t1) = (t1, t0 - q * t1);
    }
    if t0 < 0 { t0 += m as i64; }
    t0 as u64
}

// a*b mod q for q < 2^MAX_Q_BITS using 24-bit chunks of b so every intermediate fits in u128
fn mulmod_big(a: u128, b: u128, q: u128) -> u128 {
    let mut res: u128 = 0;
    let mut shift = 96i32;
    while shift >= 0 {
        let chunk = (b >> shift) & 0xFF_FFFF;
        res = (res << 24) % q;
        res = (res + a * chunk) % q;
        shift -= 24;
    }
    res
}

// 2^p mod q (left-to-right; squaring only, doubling is a shift)
fn two_pow_mod(p: u64, q: u128) -> u128 {
    let mut r: u128 = 1;
    let top = 63 - p.leading_zeros() as i32;
    for i in (0..=top).rev() {
        r = if q < (1u128 << 64) { (r * r) % q } else { mulmod_big(r, r, q) };
        if (p >> i) & 1 == 1 {
            r <<= 1;
            if r >= q { r -= q; }
        }
    }
    r
}

/// Smallest k in [k_lo, k_hi] where q = 2kp+1 divides 2^p-1, as (q, k); None if no such k.
#[pyfunction]
#[pyo3(signature = (p, k_lo, k_hi))]
fn tf_window(p: u64, k_lo: u64, k_hi: u64) -> PyResult<Option<(u128, u64)>> {
    if p < 3 || p & 1 == 0 {
        return Err(PyValueError::new_err("p must be an odd prime >= 3"));
    }
    let k_lo = k_lo.max(1);
    if k_lo > k_hi { return Ok(None); }

    let q_max: u128 = 2 * (k_hi as u128) * (p as u128) + 1;
    if q_max >= (1u128 << MAX_Q_BITS) {
        return Err(PyValueError::new_err("2*k_hi*p+1 must be < 2^100"));
    }
    let mp: Option<u128> = if p < 120 { Some((1u128 << p) - 1) } else { None };

    let width = k_hi - k_lo + 1;
    let s_limit = isqrt_u128(q_max)
        .min(SIEVE_CAP as u128)
        .min(4 * width as u128 + 100) as u64;
    let (all, n) = cached_primes(s_limit);
    let primes = &all[..n];
    let pu = p as u64;
    let two_p = 2 * pu;

    // next[i] = first k >= k_lo with s | 2kp+1 (u64::MAX = never)
    let mut next: Vec<u64> = Vec::with_capacity(primes.len());
    for &s32 in primes {
        let s = s32 as u64;
        if s == 2 || s == pu { next.push(u64::MAX); continue; }
        let k0 = s - inv_mod(two_p % s, s);
        let mut first = if k0 >= k_lo { k0 } else { k0 + (k_lo - k0).div_ceil(s) * s };
        // q == s is prime itself, not a composite marker
        if 2 * (k0 as u128) * (p as u128) + 1 == s as u128 && first == k0 { first += s; }
        next.push(first);
    }

    let mut alive = vec![1u8; BLOCK as usize];
    let mut bs = k_lo;
    while bs <= k_hi {
        let be = (bs + BLOCK - 1).min(k_hi);
        let len = (be - bs + 1) as usize;
        alive[..len].fill(1);

        for (i, &s32) in primes.iter().enumerate() {
            let s = s32 as u64;
            let mut k = next[i];
            if k > be { continue; }
            while k <= be {
                alive[(k - bs) as usize] = 0;
                k += s;
            }
            next[i] = k;
        }

        for off in 0..len {
            if alive[off] == 0 { continue; }
            let k = bs + off as u64;
            let q = 2u128 * (k as u128) * (p as u128) + 1;
            if let Some(m) = mp { if q >= m { return Ok(None); } }
            let r8 = q & 7;
            if r8 != 1 && r8 != 7 { continue; }
            if two_pow_mod(p, q) == 1 {
                return Ok(Some((q, k)));
            }
        }
        bs = be + 1;
    }
    Ok(None)
}

#[pymodule]
fn p_hunter_tf(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(tf_window, m)?)?;
    Ok(())
}
