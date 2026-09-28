//! Proven signs of the polynomials behind the separation screen.
//!
//! A breach window can be far narrower than the spacing of doubles at the instant it
//! happens. Two aircraft passing exactly through each other 1 s from now, against a
//! minimum of 1e-15 m, are inside it for 2e-17 s, and both rounded ends of that window
//! land on 1.0. Deciding from those two ends reads the window as empty. The screen
//! instead decides from the signs of polynomials in its inputs, and a sign does not
//! depend on how narrow the window is.
//!
//! Each sign is first attempted in midpoint-radius arithmetic: every operation widens a
//! radius that provably encloses the exact value, rounding the radius up. When the
//! enclosure straddles zero, the polynomial is evaluated again exactly, as a sum of
//! non-overlapping doubles: Dekker's products and Shewchuk's expansion sums. That runs on
//! a copy of the inputs scaled by powers of two, which changes no sign. Dekker's product
//! is exact only where it neither overflows nor underflows, so every product is checked
//! to lie between 2^-960 and 2^1000; if one does not, the sign is reported as unknown
//! rather than guessed.

use std::cmp::Ordering;

/// The inputs of every predicate, in this order: relative position east, north, up;
/// relative velocity east, north, up; horizontal and vertical minima; horizon.
pub type Inputs = [f64; 9];

pub const PE: usize = 0;
pub const PN: usize = 1;
pub const PU: usize = 2;
pub const VE: usize = 3;
pub const VN: usize = 4;
pub const VU: usize = 5;
pub const H: usize = 6;
pub const V: usize = 7;
pub const T: usize = 8;

const LENGTHS: [usize; 5] = [PE, PN, PU, H, V];
const SPEEDS: [usize; 3] = [VE, VN, VU];

/// Arithmetic a predicate can be written once against, and evaluated in either form.
pub trait Num: Sized {
    fn of(x: f64) -> Self;
    fn add(&self, o: &Self) -> Self;
    fn sub(&self, o: &Self) -> Self;
    fn mul(&self, o: &Self) -> Self;
    /// The sign, if it is proven.
    fn sign(&self) -> Option<Ordering>;
}

/// The sign of the polynomial `f` at `x`, proven, or `None`. `uses` lists the inputs `f`
/// reads; the exact evaluation scales those, and only those.
pub fn sign(
    x: &Inputs,
    uses: &[usize],
    f_ball: fn(&[Ball; 9]) -> Ball,
    f_exact: fn(&[Exact; 9]) -> Exact,
) -> Option<Ordering> {
    if let Some(s) = f_ball(&x.map(Ball::of)).sign() {
        return Some(s);
    }
    let x = scaled(x, uses).unwrap_or(*x);
    f_exact(&x.map(Exact::of)).sign()
}

/// 2^-52, which bounds the relative error of one rounding in terms of its result.
const U: f64 = f64::EPSILON;
/// The smallest subnormal, which bounds what a product can lose by underflowing.
const ETA: f64 = f64::from_bits(1);

fn up(x: f64) -> f64 {
    x.next_up()
}

/// A value known to lie within `r` of `m`. Each operation computes the radius with every
/// rounding taken upward, so the enclosure holds whatever the rounding of `m` did.
#[derive(Clone, Copy, Debug)]
pub struct Ball {
    m: f64,
    r: f64,
}

impl Num for Ball {
    fn of(x: f64) -> Self {
        Ball { m: x, r: 0.0 }
    }

    fn add(&self, o: &Self) -> Self {
        // A sum never underflows inexactly, so only the relative term applies.
        let m = self.m + o.m;
        let r = up(up(self.r + o.r) + up(U * m.abs()));
        Ball { m, r }
    }

    fn sub(&self, o: &Self) -> Self {
        let m = self.m - o.m;
        let r = up(up(self.r + o.r) + up(U * m.abs()));
        Ball { m, r }
    }

    fn mul(&self, o: &Self) -> Self {
        // (a + da)(b + db) − ab = a·db + b·da + da·db, then the rounding of ab itself,
        // which can underflow by up to half the smallest subnormal.
        let m = self.m * o.m;
        let spread = up(up(up(self.m.abs() * o.r) + up(o.m.abs() * self.r)) + up(self.r * o.r));
        let r = up(up(spread + up(U * m.abs())) + ETA);
        Ball { m, r }
    }

    fn sign(&self) -> Option<Ordering> {
        if !(self.m.is_finite() && self.r.is_finite()) {
            None
        } else if self.m > self.r {
            Some(Ordering::Greater)
        } else if -self.m > self.r {
            Some(Ordering::Less)
        } else {
            None
        }
    }
}

/// A number held exactly as a sum of non-overlapping doubles in increasing magnitude,
/// zeros dropped, or marked unusable once any step could not be exact.
#[derive(Clone, Debug)]
pub struct Exact {
    c: Vec<f64>,
    ok: bool,
}

fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let x = a + b;
    let bv = x - a;
    let av = x - bv;
    (x, (a - av) + (b - bv))
}

/// 2^995: above this, Dekker's splitting constant times the operand overflows.
const SPLIT_MAX: f64 = f64::from_bits((1023 + 995) << 52);
/// 2^1000 and 2^-960: between these, a product and its error are both exactly
/// representable. The lower bound keeps the operands' exponents summing to at least
/// −970, where the smallest partial product still has its last bit above 2^-1074.
const PRODUCT_MAX: f64 = f64::from_bits((1023 + 1000) << 52);
const PRODUCT_MIN: f64 = f64::from_bits((1023 - 960) << 52);

fn split(a: f64) -> (f64, f64) {
    let c = 134_217_729.0 * a;
    let hi = c - (c - a);
    (hi, a - hi)
}

/// `a·b` as a rounded product and its exact error, or `None` outside the exact range.
fn two_prod(a: f64, b: f64) -> Option<(f64, f64)> {
    let x = a * b;
    if a == 0.0 || b == 0.0 {
        return Some((x, 0.0));
    }
    let fits = a.abs() <= SPLIT_MAX
        && b.abs() <= SPLIT_MAX
        && x.abs() <= PRODUCT_MAX
        && x.abs() >= PRODUCT_MIN;
    if !fits {
        return None;
    }
    let (ah, al) = split(a);
    let (bh, bl) = split(b);
    let err = ((x - ah * bh) - al * bh) - ah * bl;
    Some((x, al * bl - err))
}

impl Exact {
    fn finish(c: Vec<f64>, ok: bool) -> Exact {
        let ok = ok && c.iter().all(|v| v.is_finite());
        Exact { c, ok }
    }

    /// This expansion plus one double.
    fn grow(&self, b: f64) -> Exact {
        let mut out = Vec::with_capacity(self.c.len() + 1);
        let mut q = b;
        for &e in &self.c {
            let (s, h) = two_sum(q, e);
            if h != 0.0 {
                out.push(h);
            }
            q = s;
        }
        if q != 0.0 {
            out.push(q);
        }
        Exact::finish(out, self.ok)
    }

    fn sum(&self, o: &Exact, negate: bool) -> Exact {
        let mut acc = Exact {
            c: self.c.clone(),
            ok: self.ok && o.ok,
        };
        for &f in &o.c {
            acc = acc.grow(if negate { -f } else { f });
        }
        acc
    }

    /// This expansion times one double.
    fn scale(&self, b: f64) -> Exact {
        let mut out = Vec::with_capacity(2 * self.c.len());
        let mut parts = self.c.iter();
        let Some(&first) = parts.next() else {
            return Exact::finish(out, self.ok);
        };
        let Some((mut q, h)) = two_prod(first, b) else {
            return Exact::finish(out, false);
        };
        if h != 0.0 {
            out.push(h);
        }
        for &e in parts {
            let Some((hi, lo)) = two_prod(e, b) else {
                return Exact::finish(out, false);
            };
            let (s, h) = two_sum(q, lo);
            if h != 0.0 {
                out.push(h);
            }
            let (s, h) = two_sum(hi, s);
            if h != 0.0 {
                out.push(h);
            }
            q = s;
        }
        if q != 0.0 {
            out.push(q);
        }
        Exact::finish(out, self.ok)
    }
}

impl Num for Exact {
    fn of(x: f64) -> Self {
        let c = if x == 0.0 { Vec::new() } else { vec![x] };
        Exact::finish(c, true)
    }

    fn add(&self, o: &Self) -> Self {
        self.sum(o, false)
    }

    fn sub(&self, o: &Self) -> Self {
        self.sum(o, true)
    }

    fn mul(&self, o: &Self) -> Self {
        let mut acc = Exact {
            c: Vec::new(),
            ok: self.ok && o.ok,
        };
        for &f in &o.c {
            acc = acc.sum(&self.scale(f), false);
        }
        acc
    }

    fn sign(&self) -> Option<Ordering> {
        // The largest component carries the sign of the whole sum.
        if !self.ok {
            return None;
        }
        Some(match self.c.last() {
            None => Ordering::Equal,
            Some(&v) if v > 0.0 => Ordering::Greater,
            Some(_) => Ordering::Less,
        })
    }
}

/// The binary exponent of a nonzero finite double.
fn exponent(x: f64) -> i32 {
    let bits = x.to_bits() & !(1u64 << 63);
    let biased = (bits >> 52) as i32;
    if biased == 0 {
        -1074 + (63 - bits.leading_zeros() as i32)
    } else {
        biased - 1023
    }
}

/// `x·2^k`, if that is exact.
fn times_pow2(x: f64, k: i32) -> Option<f64> {
    let apply = |mut y: f64, mut k: i32| {
        while k != 0 {
            let step = k.clamp(-1000, 1000);
            y *= f64::from_bits(((1023 + step) as u64) << 52);
            k -= step;
        }
        y
    };
    let y = apply(x, k);
    (y.is_finite() && apply(y, -k) == x).then_some(y)
}

/// The inputs `uses` names, with lengths scaled by one power of two and velocities by
/// another, and the horizon by their ratio, so that each unit's magnitudes are centred
/// on 1. Every predicate is homogeneous in each unit, so no sign changes. `None` if a
/// scaled value would not be exact.
fn scaled(x: &Inputs, uses: &[usize]) -> Option<Inputs> {
    let centre = |class: &[usize]| {
        let (lo, hi) = class
            .iter()
            .filter(|&&i| uses.contains(&i) && x[i] != 0.0)
            .map(|&i| exponent(x[i]))
            .fold((i32::MAX, i32::MIN), |(lo, hi), e| (lo.min(e), hi.max(e)));
        if lo > hi {
            0
        } else {
            -(lo + hi) / 2
        }
    };
    let (a, c) = (centre(&LENGTHS), centre(&SPEEDS));
    let mut out = *x;
    for &i in uses {
        let k = if LENGTHS.contains(&i) {
            a
        } else if SPEEDS.contains(&i) {
            c
        } else {
            a - c
        };
        out[i] = times_pow2(x[i], k)?;
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// (pe·vn − pn·ve)² − H²(ve² + vn²), the degree-4 test for a track passing inside a
    /// minimum.
    fn cross_test<N: Num>(x: &[N; 9]) -> N {
        let cross = x[PE].mul(&x[VN]).sub(&x[PN].mul(&x[VE]));
        let speed2 = x[VE].mul(&x[VE]).add(&x[VN].mul(&x[VN]));
        cross.mul(&cross).sub(&x[H].mul(&x[H]).mul(&speed2))
    }

    fn truth(p: [i64; 5]) -> Ordering {
        let [pe, pn, ve, vn, h] = p.map(i128::from);
        let cross = pe * vn - pn * ve;
        (cross * cross).cmp(&(h * h * (ve * ve + vn * vn)))
    }

    struct Lcg(u64);
    impl Lcg {
        fn below(&mut self, n: u64) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 33) % n
        }
        fn int(&mut self, bits: u32) -> i64 {
            let m = 1i64 << bits;
            self.below(2 * m as u64) as i64 - m
        }
    }

    #[test]
    fn signs_match_exact_integer_arithmetic_at_every_scale() {
        // Random integers up to 2^25, and tangencies from the 3-4-5 triangle: with
        // v = (3m, 4m) and p = (3τ + 4d, 4τ − 3d), the track passes 5|d| from the origin,
        // and H = 5|d| is tangent. With m, τ and d up to 2^26 the products need up to 57
        // bits, so both the rounded products and their exact errors decide the answer. The
        // same cases are then scaled by powers of two up to 2^±900, which must change no
        // sign.
        let mut rng = Lcg(0x5e9a_7a7e);
        let (mut exact_needed, mut unscaled) = (0, 0);
        for case in 0..20_000 {
            let p = if case % 2 == 0 {
                [
                    rng.int(25),
                    rng.int(25),
                    rng.int(25),
                    rng.int(25),
                    rng.int(25).abs(),
                ]
            } else {
                let (m, tau, d) = (rng.int(26), rng.int(26), rng.int(26));
                let h = (5 * d).abs() + [-1, 0, 0, 1][rng.below(4) as usize];
                [3 * tau + 4 * d, 4 * tau - 3 * d, 3 * m, 4 * m, h]
            };
            let want = truth(p);
            let base: Inputs = {
                let mut x = [0.0; 9];
                (x[PE], x[PN], x[VE], x[VN], x[H]) = (
                    p[0] as f64,
                    p[1] as f64,
                    p[2] as f64,
                    p[3] as f64,
                    p[4] as f64,
                );
                x
            };
            let uses = [PE, PN, VE, VN, H];
            if let Some(s) = cross_test(&base.map(Ball::of)).sign() {
                assert_eq!(s, want, "the enclosure claimed a wrong sign: {p:?}");
            } else {
                exact_needed += 1;
            }
            let got = sign(&base, &uses, cross_test::<Ball>, cross_test::<Exact>);
            assert_eq!(got, Some(want), "{p:?}");
            if cross_test(&base.map(Exact::of)).sign() == Some(want) {
                unscaled += 1;
            }
            let kl = (rng.int(10) as i32).clamp(-900, 900);
            let kv = (rng.int(10) as i32).clamp(-900, 900);
            let mut x = base;
            for &i in &[PE, PN, H] {
                x[i] = times_pow2(x[i], kl).unwrap();
            }
            for &i in &[VE, VN] {
                x[i] = times_pow2(x[i], kv).unwrap();
            }
            let got = sign(&x, &uses, cross_test::<Ball>, cross_test::<Exact>);
            assert_eq!(got, Some(want), "{p:?} scaled by 2^{kl}, 2^{kv}");
        }
        // Neither path may be vacuous.
        assert!(
            exact_needed > 1_000,
            "{exact_needed} cases needed the exact path"
        );
        assert_eq!(unscaled, 20_000);
    }

    #[test]
    fn a_product_outside_the_exact_range_is_unknown_not_guessed() {
        let tiny = f64::from_bits((1023 - 500) << 52);
        let x = Exact::of(tiny).mul(&Exact::of(tiny));
        assert_eq!(x.sign(), None);
        let big = f64::from_bits((1023 + 600) << 52);
        assert_eq!(Exact::of(big).mul(&Exact::of(big)).sign(), None);
        assert_eq!(
            Exact::of(0.0).mul(&Exact::of(big)).sign(),
            Some(Ordering::Equal)
        );
    }
}
