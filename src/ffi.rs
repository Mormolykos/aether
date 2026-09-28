//! C ABI surface, so the estimator and the conjunction screen can be called from an
//! existing C or C++ system without rewriting either side.
//!
//! Rules held on this boundary:
//!   * `#[repr(C)]` on everything that crosses it, so the layout is the C layout.
//!   * every pointer is null-checked before it is read; a bad pointer returns an error
//!     code rather than dereferencing.
//!   * no panic is allowed to unwind into the caller — unwinding across an FFI boundary
//!     is undefined behaviour, so every entry point is panic-free by construction.
//!   * no allocation crosses the boundary, so there is no question of which allocator
//!     frees what.
//!   * every input is checked against the declared numerical domain, and every result
//!     is checked before it is written: success means four finite numbers.
//!
//! Matching C header:
//!
//! ```c
//! typedef struct { double e, n, u, ve, vn, vu; } aether_state_t;
//! typedef struct { double t_cpa, horiz_m, vert_m, closing_ms; } aether_cpa_t;
//!
//! #define AETHER_OK         0
//! #define AETHER_ERR_NULL  -1  /* a pointer was null */
//! #define AETHER_ERR_NAN   -2  /* an input or the horizon was NaN, infinite or negative */
//! #define AETHER_ERR_RANGE -3  /* a finite input exceeded 1e9 in magnitude, or the
//!                                 result could not be represented */
//!
//! int aether_cpa(const aether_state_t *a, const aether_state_t *b,
//!                double horizon_s, aether_cpa_t *out);
//! ```

use crate::conjunction::closest_approach_s;
use crate::domain::{within, ENVELOPE};
use crate::geo::Enu;
use std::os::raw::c_int;

pub const AETHER_OK: c_int = 0;
pub const AETHER_ERR_NULL: c_int = -1;
pub const AETHER_ERR_NAN: c_int = -2;
/// A finite input outside the declared envelope (|x| > 1e9, SI units), or a result that
/// could not be represented. Distinct from `AETHER_ERR_NAN` so a caller can tell a
/// corrupted value from an out-of-range one.
pub const AETHER_ERR_RANGE: c_int = -3;

/// Position and velocity in a local ENU frame, metres and metres per second.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct AetherState {
    pub e: f64,
    pub n: f64,
    pub u: f64,
    pub ve: f64,
    pub vn: f64,
    pub vu: f64,
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct AetherCpa {
    pub t_cpa: f64,
    pub horiz_m: f64,
    pub vert_m: f64,
    pub closing_ms: f64,
}

/// Closest point of approach between two constant-velocity states.
///
/// A closest-approach primitive, not a separation verdict. `horiz_m` and `vert_m` are the
/// separations at the instant of least 3-D distance; a pair can be inside both minima at
/// another instant while outside one of them at this one. Deciding loss of separation
/// from these fields alone repeats that error. `conjunction::pair_cpa` decides it from
/// the violation intervals instead.
///
/// Returns `AETHER_OK` and fills `out` on success. `out` is untouched on any error.
///
/// # Safety
/// `a`, `b` and `out` must each be either null or a valid, aligned pointer to a single
/// initialised value of the corresponding type. Null is handled and reported.
#[no_mangle]
pub unsafe extern "C" fn aether_cpa(
    a: *const AetherState,
    b: *const AetherState,
    horizon_s: f64,
    out: *mut AetherCpa,
) -> c_int {
    if a.is_null() || b.is_null() || out.is_null() {
        return AETHER_ERR_NULL;
    }
    let (a, b) = (&*a, &*b);
    if !horizon_s.is_finite() || horizon_s < 0.0 {
        return AETHER_ERR_NAN;
    }

    match cpa(a, b, horizon_s) {
        Ok(result) => {
            out.write(result);
            AETHER_OK
        }
        Err(code) => code,
    }
}

/// The same closest approach that `conjunction::pair_cpa` reports, through the same
/// function, on bare states. Kept as safe Rust so it can be unit tested without any
/// unsafe block in the test.
///
/// Checking only that the inputs are finite is not enough: two finite states of order
/// 1e200 overflow the squared velocity, and the time of closest approach becomes NaN. So
/// every input is held to the envelope first, which makes overflow unreachable, and the
/// result is checked as well, so that success can only ever carry finite numbers.
fn cpa(a: &AetherState, b: &AetherState, horizon_s: f64) -> Result<AetherCpa, c_int> {
    let inputs = [
        a.e, a.n, a.u, a.ve, a.vn, a.vu, b.e, b.n, b.u, b.ve, b.vn, b.vu,
    ];
    if inputs.iter().any(|v| !v.is_finite()) {
        return Err(AETHER_ERR_NAN);
    }
    if !inputs.iter().all(|&v| within(v, -ENVELOPE, ENVELOPE)) {
        return Err(AETHER_ERR_RANGE);
    }

    let dp = Enu {
        e: a.e - b.e,
        n: a.n - b.n,
        u: a.u - b.u,
    };
    let dv = Enu {
        e: a.ve - b.ve,
        n: a.vn - b.vn,
        u: a.vu - b.vu,
    };
    let t = closest_approach_s(dp, dv, horizon_s);
    let at = dp + dv * t;
    let result = AetherCpa {
        t_cpa: t,
        horiz_m: at.horiz(),
        vert_m: at.u.abs(),
        closing_ms: dv.e.hypot(dv.n).hypot(dv.u),
    };
    let fields = [
        result.t_cpa,
        result.horiz_m,
        result.vert_m,
        result.closing_ms,
    ];
    if fields.iter().all(|v| v.is_finite()) {
        Ok(result)
    } else {
        Err(AETHER_ERR_RANGE)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn state(e: f64, n: f64, u: f64, ve: f64, vn: f64, vu: f64) -> AetherState {
        AetherState {
            e,
            n,
            u,
            ve,
            vn,
            vu,
        }
    }

    #[test]
    fn head_on_pair_closes_at_the_midpoint() {
        // 20 km apart on the north axis, 100 m/s each toward the other.
        let a = state(0.0, 10_000.0, 0.0, 0.0, -100.0, 0.0);
        let b = state(0.0, -10_000.0, 0.0, 0.0, 100.0, 0.0);
        let r = cpa(&a, &b, 300.0).unwrap();
        assert!((r.t_cpa - 100.0).abs() < 1e-6, "t_cpa {}", r.t_cpa);
        assert!(
            r.horiz_m < 1e-6,
            "should pass through zero, got {}",
            r.horiz_m
        );
        assert!((r.closing_ms - 200.0).abs() < 1e-6);
    }

    #[test]
    fn a_past_approach_is_clamped_to_now() {
        // Already separating: the mathematical minimum is in the past.
        let a = state(0.0, 100.0, 0.0, 0.0, 100.0, 0.0);
        let b = state(0.0, -100.0, 0.0, 0.0, -100.0, 0.0);
        let r = cpa(&a, &b, 300.0).unwrap();
        assert_eq!(r.t_cpa, 0.0);
        assert!((r.horiz_m - 200.0).abs() < 1e-6);
    }

    #[test]
    fn co_velocity_pair_never_converges() {
        let a = state(0.0, 0.0, 0.0, 250.0, 0.0, 0.0);
        let b = state(5_000.0, 0.0, 0.0, 250.0, 0.0, 0.0);
        let r = cpa(&a, &b, 300.0).unwrap();
        assert_eq!(r.t_cpa, 0.0);
        assert!((r.horiz_m - 5_000.0).abs() < 1e-6);
        assert_eq!(r.closing_ms, 0.0);
    }

    #[test]
    fn vertical_separation_is_reported_separately() {
        let a = state(0.0, 10_000.0, 3_000.0, 0.0, -100.0, 0.0);
        let b = state(0.0, -10_000.0, 3_400.0, 0.0, 100.0, 0.0);
        let r = cpa(&a, &b, 300.0).unwrap();
        assert!(r.horiz_m < 1e-6);
        assert!((r.vert_m - 400.0).abs() < 1e-6);
    }

    #[test]
    fn null_pointers_are_reported_not_dereferenced() {
        let s = state(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let mut out = AetherCpa::default();
        unsafe {
            assert_eq!(
                aether_cpa(std::ptr::null(), &s, 60.0, &mut out),
                AETHER_ERR_NULL
            );
            assert_eq!(
                aether_cpa(&s, std::ptr::null(), 60.0, &mut out),
                AETHER_ERR_NULL
            );
            assert_eq!(
                aether_cpa(&s, &s, 60.0, std::ptr::null_mut()),
                AETHER_ERR_NULL
            );
        }
    }

    #[test]
    fn nan_input_is_rejected_rather_than_propagated() {
        let good = state(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let bad = state(f64::NAN, 0.0, 0.0, 0.0, 0.0, 0.0);
        let mut out = AetherCpa::default();
        unsafe {
            assert_eq!(aether_cpa(&good, &bad, 60.0, &mut out), AETHER_ERR_NAN);
            assert_eq!(aether_cpa(&good, &good, f64::NAN, &mut out), AETHER_ERR_NAN);
        }
        // out must be untouched on the error path.
        assert_eq!(out.t_cpa, 0.0);
    }

    /// Output pre-filled with a bit pattern no computation produces, so "untouched" can
    /// be checked exactly rather than by value.
    fn sentinel() -> AetherCpa {
        let s = f64::from_bits(0x7ff4_dead_beef_0001);
        AetherCpa {
            t_cpa: s,
            horiz_m: s,
            vert_m: s,
            closing_ms: s,
        }
    }

    fn untouched(out: &AetherCpa) -> bool {
        let s = sentinel();
        [out.t_cpa, out.horiz_m, out.vert_m, out.closing_ms]
            .iter()
            .zip([s.t_cpa, s.horiz_m, s.vert_m, s.closing_ms])
            .all(|(a, b)| a.to_bits() == b.to_bits())
    }

    #[test]
    fn finite_inputs_that_overflow_are_refused_and_leave_output_untouched() {
        // Audit regression. Every component is finite, so the input check passed; the
        // squared velocity and the dot product overflow, time becomes NaN, and the old
        // entry point wrote that and returned success.
        let a = state(1e200, 0.0, 0.0, -1e200, 0.0, 0.0);
        let b = state(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let mut out = sentinel();
        let rc = unsafe { aether_cpa(&a, &b, 300.0, &mut out) };
        assert_ne!(
            rc, AETHER_OK,
            "overflowing arithmetic reported success: {out:?}"
        );
        assert_eq!(rc, AETHER_ERR_RANGE);
        assert!(
            untouched(&out),
            "an error path wrote to the caller: {out:?}"
        );
    }

    #[test]
    fn a_very_slow_approach_is_not_reported_as_already_closest() {
        // Audit regression. 1 m apart, closing at 1e-5 m/s, over a 100,000 s horizon: the
        // closest approach is at 100,000 s and 0 m. A threshold that treated any squared
        // speed under 1e-9 as parallel reported t = 0 and 1 m instead.
        let a = state(1.0, 0.0, 0.0, -1e-5, 0.0, 0.0);
        let b = state(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let r = cpa(&a, &b, 100_000.0).unwrap();
        assert!((r.t_cpa - 100_000.0).abs() < 1e-6, "{r:?}");
        assert!(r.horiz_m < 1e-9, "{r:?}");

        // Neighbouring scales: speeds from 1e-3 down to 1e-12 m/s, each over the horizon
        // it takes to close 1 m, against a direct evaluation of the same instant.
        for k in 3..=12 {
            let speed = 10f64.powi(-k);
            let horizon = 1.0 / speed;
            let a = state(1.0, 0.0, 0.0, -speed, 0.0, 0.0);
            let r = cpa(&a, &b, horizon).unwrap();
            assert!(
                (r.t_cpa / horizon - 1.0).abs() < 1e-9 && r.horiz_m < 1e-6,
                "speed {speed}: {r:?}"
            );
            assert!(
                (r.closing_ms / speed - 1.0).abs() < 1e-12,
                "speed {speed}: {r:?}"
            );
        }
        // One geometry at every scale the envelope admits: s metres apart, closing at s
        // metres per second, meets after 1 s. Below about 1e-154 the squares of s
        // underflow, and a quotient formed from them is 0/0.
        for s in [1e-100, 1e-160, 1e-200, 1e-300, f64::MIN_POSITIVE, 5e-324] {
            let r = cpa(&state(s, 0.0, 0.0, -s, 0.0, 0.0), &b, 10.0).unwrap();
            assert!(
                r.t_cpa == 1.0 && r.horiz_m == 0.0 && r.closing_ms == s,
                "scale {s}: {r:?}"
            );
        }
        // Exactly equal velocities have no closest approach later than now.
        let still = cpa(&b, &state(3.0, 4.0, 0.0, 0.0, 0.0, 0.0), 1e9).unwrap();
        assert_eq!((still.t_cpa, still.horiz_m), (0.0, 5.0));
    }

    #[test]
    fn inputs_at_the_edge_of_the_envelope_still_succeed() {
        // The envelope refuses what cannot be real, not what is merely large: two states
        // a billion metres apart, closing at a billion metres per second, are computed.
        let a = state(1e9, -1e9, 1e9, -1e9, 1e9, -1e9);
        let b = state(-1e9, 1e9, -1e9, 1e9, -1e9, 1e9);
        let mut out = sentinel();
        let rc = unsafe { aether_cpa(&a, &b, 300.0, &mut out) };
        assert_eq!(rc, AETHER_OK);
        assert!(
            out.t_cpa.is_finite() && out.closing_ms.is_finite(),
            "{out:?}"
        );
        let over = state(1.000_000_1e9, 0.0, 0.0, 0.0, 0.0, 0.0);
        assert_eq!(
            unsafe { aether_cpa(&over, &b, 300.0, &mut out) },
            AETHER_ERR_RANGE
        );
    }

    #[test]
    fn every_call_either_succeeds_finitely_or_fails_without_writing() {
        // The property across magnitudes: ordinary, large, and at the edge of overflow,
        // in every field and with every sign, over ordinary and absurd horizons.
        let mags = [
            0.0, 1.0, 250.0, 1e7, 1e9, 1e12, 1e100, 1e154, 1e155, 1e200, 1e308,
        ];
        let horizons = [0.0, 300.0, 1e9, 1e300, f64::MAX];
        let mut calls = 0;
        for (i, &m) in mags.iter().enumerate() {
            for field in 0..6 {
                for sign in [1.0, -1.0] {
                    for &h in &horizons {
                        let mut v = [100.0, -2_000.0, 30.0, -150.0, 80.0, 5.0];
                        v[field] = sign * m;
                        // A second large field, so products of two large numbers occur.
                        v[(field + 3) % 6] = -sign * mags[(i + 3) % mags.len()];
                        let a = state(v[0], v[1], v[2], v[3], v[4], v[5]);
                        let b = state(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                        let mut out = sentinel();
                        let rc = unsafe { aether_cpa(&a, &b, h, &mut out) };
                        calls += 1;
                        if rc == AETHER_OK {
                            let all = [out.t_cpa, out.horiz_m, out.vert_m, out.closing_ms];
                            assert!(
                                all.iter().all(|x| x.is_finite()),
                                "success with a non-finite result for {v:?}, horizon {h}: {out:?}"
                            );
                        } else {
                            assert!(untouched(&out), "error {rc} wrote for {v:?}: {out:?}");
                        }
                    }
                }
            }
        }
        assert!(calls > 600);
    }

    #[test]
    fn ffi_call_succeeds_through_the_real_entry_point() {
        let a = state(0.0, 10_000.0, 0.0, 0.0, -100.0, 0.0);
        let b = state(0.0, -10_000.0, 0.0, 0.0, 100.0, 0.0);
        let mut out = AetherCpa::default();
        let rc = unsafe { aether_cpa(&a, &b, 300.0, &mut out) };
        assert_eq!(rc, AETHER_OK);
        assert!((out.t_cpa - 100.0).abs() < 1e-6);
    }
}
