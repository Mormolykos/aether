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
//!
//! Matching C header:
//!
//! ```c
//! typedef struct { double e, n, u, ve, vn, vu; } aether_state_t;
//! typedef struct { double t_cpa, horiz_m, vert_m, closing_ms; } aether_cpa_t;
//!
//! int aether_cpa(const aether_state_t *a, const aether_state_t *b,
//!                double horizon_s, aether_cpa_t *out);
//! ```

use std::os::raw::c_int;

pub const AETHER_OK: c_int = 0;
pub const AETHER_ERR_NULL: c_int = -1;
pub const AETHER_ERR_NAN: c_int = -2;

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

    let result = match cpa(a, b, horizon_s) {
        Some(r) => r,
        None => return AETHER_ERR_NAN,
    };
    out.write(result);
    AETHER_OK
}

/// The same closed form as `conjunction::pair_cpa`, on bare states. Kept as safe Rust
/// so it can be unit tested without any unsafe block in the test.
fn cpa(a: &AetherState, b: &AetherState, horizon_s: f64) -> Option<AetherCpa> {
    let (dpe, dpn, dpu) = (a.e - b.e, a.n - b.n, a.u - b.u);
    let (dve, dvn, dvu) = (a.ve - b.ve, a.vn - b.vn, a.vu - b.vu);

    let fields = [dpe, dpn, dpu, dve, dvn, dvu];
    if fields.iter().any(|v| !v.is_finite()) {
        return None;
    }

    let dvv = dve * dve + dvn * dvn + dvu * dvu;
    let t = if dvv < 1e-9 {
        0.0
    } else {
        (-(dpe * dve + dpn * dvn + dpu * dvu) / dvv).clamp(0.0, horizon_s)
    };

    let (e, n, u) = (dpe + dve * t, dpn + dvn * t, dpu + dvu * t);
    Some(AetherCpa {
        t_cpa: t,
        horiz_m: e.hypot(n),
        vert_m: u.abs(),
        closing_ms: (dve * dve + dvn * dvn + dvu * dvu).sqrt(),
    })
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
