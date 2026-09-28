//! The declared numerical domain.
//!
//! The predicates for the values Aether accepts at three boundaries: a contact entering
//! `TrackStore::ingest`, together with the frame it is transformed through; a
//! configuration file; and a state handed across the C ABI. A value that fails its
//! predicate is treated as absent: refused, or for an optional field dropped. It is
//! never repaired into something plausible. Arithmetic after those boundaries is guarded
//! separately — a filter update is kept only if its result passes `track::usable` —
//! because bounded inputs do not by themselves bound every intermediate: a small
//! measurement variance and a tiny time step produce enormous, finite gains.
//!
//! # Per-quantity bounds
//!
//! | quantity | unit | accepted | basis |
//! |---|---|---|---|
//! | latitude, longitude | degree | [−90, 90], [−180, 180] | definition |
//! | altitude (contact, observer) | m | within ±1e9 | policy envelope |
//! | reported age | s | [0, 1e9] | ≥ 0 by definition; upper bound policy |
//! | ground speed, climb (optional) | m/s | within ±1e9, speed ≥ 0 | policy envelope; outside it the field is absent |
//! | track over ground (optional) | degree | [0, 360] | definition; outside it the field is absent |
//! | estimated velocity, per axis | m/s | within ±1e9 | policy envelope, checked per axis rather than as a speed |
//! | process noise q | m²/s³ | 0, or a normal double up to 1e9 | zero is a model choice; see below for the floor |
//! | measurement variance r | m² | a normal double up to 1e9 | must be > 0, see below |
//! | gate | σ | a normal double up to 1e9 | must be > 0, see below |
//! | separation minima | m | a normal double up to 1e9 | must be > 0, see below |
//! | horizon | s | 0, or a normal double up to 1e9 | zero means "now only" |
//! | C ABI position, velocity components | m, m/s | within ±1e9 | policy envelope |
//!
//! **The 1e9 envelope is a policy**, not a physical law and not a derived necessity. It is
//! set far above the regional air picture Aether is built for: aircraft positions within
//! a few hundred kilometres and speeds below 1 km/s. It is the same number in every unit
//! only for simplicity. A measurement variance above 1e9 m² (a sigma above 31.6 km) is not
//! physically impossible; this demonstrator declines it.
//!
//! **The floor at the smallest normal double (about 2.2e-308) is derived** from the number
//! format. Below it a double carries fewer than 53 significant bits — at 5e-324, one — so
//! a configured value like `3e-324` is read as 5e-324, and a variance that small makes the
//! covariance itself unrepresentable even though its square-root factor is not. Above the
//! floor, no lower bound is imposed on these quantities: the filter update and the
//! horizontal interval are formed without the squares and subtractions that failed at
//! small values (see `kalman` and `conjunction`), and are tested down to the floor.
//! Whether a very small minimum or gate is operationally sensible is the operator's
//! policy, not a numerical question.

use crate::geo::{Geodetic, FPM_TO_MS, KT_TO_MS};
use crate::ingest::Contact;

/// Largest magnitude accepted for the quantities in the module table, in SI units. A
/// policy; see the module documentation.
pub const ENVELOPE: f64 = 1.0e9;

/// `x` is a number in `[lo, hi]`. NaN fails every comparison, so it fails this, and an
/// infinity fails it whenever the bounds are finite.
pub fn within(x: f64, lo: f64, hi: f64) -> bool {
    x >= lo && x <= hi
}

/// `x` is a positive normal double no larger than `hi`: strictly positive, and carrying
/// its full 53 significant bits.
pub fn positive(x: f64, hi: f64) -> bool {
    within(x, f64::MIN_POSITIVE, hi)
}

/// Why a geodetic position is unusable, or `None` if it is usable.
pub fn geodetic(g: Geodetic) -> Option<&'static str> {
    if !within(g.lat_deg, -90.0, 90.0) {
        Some("latitude must be a number from -90 to 90 degrees")
    } else if !within(g.lon_deg, -180.0, 180.0) {
        Some("longitude must be a number from -180 to 180 degrees")
    } else if !within(g.alt_m, -ENVELOPE, ENVELOPE) {
        Some("altitude must be a number of metres within ±1e9")
    } else {
        None
    }
}

/// Why a contact cannot enter the tracker, or `None` if it can.
///
/// Only the fields a track cannot exist without are checked here. The optional kinematic
/// fields go through [`ground_speed_ms`], [`bearing_rad`] and [`climb_ms`], which turn an
/// unusable value into an absent one: a report with a nonsense ground speed still has a
/// perfectly good position, and absent is a state the tracker already handles.
pub fn contact(c: &Contact) -> Option<&'static str> {
    geodetic(c.geo).or_else(|| {
        // A report cannot be observed after it arrived, so its age is never negative.
        (!within(c.age_s, 0.0, ENVELOPE)).then_some("age must be a number of seconds from 0 to 1e9")
    })
}

/// A reported ground speed in knots, as metres per second, if present and usable.
pub fn ground_speed_ms(kt: Option<f64>) -> Option<f64> {
    kt.map(|k| k * KT_TO_MS)
        .filter(|&v| within(v, 0.0, ENVELOPE))
}

/// A reported track over the ground in degrees, as radians, if present and usable.
pub fn bearing_rad(deg: Option<f64>) -> Option<f64> {
    deg.filter(|&d| within(d, 0.0, 360.0)).map(f64::to_radians)
}

/// A reported vertical rate in feet per minute, as metres per second, if present and
/// usable.
pub fn climb_ms(fpm: Option<f64>) -> Option<f64> {
    fpm.map(|f| f * FPM_TO_MS)
        .filter(|&v| within(v, -ENVELOPE, ENVELOPE))
}

/// Why a process-noise density is unusable, or `None`. Zero is legitimate: it asserts a
/// target that never accelerates.
pub fn process_noise(q: f64) -> Option<&'static str> {
    (!(q == 0.0 || positive(q, ENVELOPE)))
        .then_some("must be 0, or a number from 2.2250738585072014e-308 to 1e9")
}

/// Why a measurement variance is unusable, or `None`. Zero is not usable. With no process
/// noise it makes the whole covariance exactly zero at the first update and the
/// innovation variance zero at the second, after which every plot is refused. With
/// process noise it runs, but it asserts a sensor with no noise at all, which is a
/// modelling choice this tracker does not make.
pub fn meas_var(r: f64) -> Option<&'static str> {
    (!positive(r, ENVELOPE)).then_some("must be a number from 2.2250738585072014e-308 to 1e9")
}

/// Why a gate width is unusable, or `None`. NaN would disable the gate outright, and zero
/// would refuse every plot a firm track is ever offered. How wide it should be is a
/// detection policy: a gate below about 3σ refuses a visible share of consistent plots,
/// and that shows in the gated count rather than being refused here.
pub fn gate_sigma(g: f64) -> Option<&'static str> {
    (!positive(g, ENVELOPE)).then_some("must be a number from 2.2250738585072014e-308 to 1e9")
}

/// Why a separation minimum is unusable, or `None`. Any positive minimum defines a
/// non-empty region and the screen's arithmetic handles it without squaring it; how small
/// a minimum is sensible is airspace policy.
pub fn separation_minimum(m: f64) -> Option<&'static str> {
    (!positive(m, ENVELOPE)).then_some("must be a number from 2.2250738585072014e-308 to 1e9")
}

/// Why a look-ahead horizon is unusable, or `None`. Zero is legitimate: screen the
/// present only.
pub fn horizon(h: f64) -> Option<&'static str> {
    (!(h == 0.0 || positive(h, ENVELOPE)))
        .then_some("must be 0, or a number of seconds from 2.2250738585072014e-308 to 1e9")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nan_and_infinity_are_outside_every_range() {
        for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(!within(x, -ENVELOPE, ENVELOPE), "{x}");
        }
        assert!(within(-0.0, 0.0, 1.0));
        assert!(within(ENVELOPE, -ENVELOPE, ENVELOPE));
    }

    #[test]
    fn an_unusable_optional_field_becomes_absent() {
        assert_eq!(ground_speed_ms(Some(1e308)), None);
        assert_eq!(ground_speed_ms(Some(-1.0)), None);
        assert_eq!(bearing_rad(Some(f64::NAN)), None);
        assert_eq!(bearing_rad(Some(361.0)), None);
        assert_eq!(climb_ms(Some(f64::INFINITY)), None);
        assert!((ground_speed_ms(Some(450.0)).unwrap() - 231.5).abs() < 0.1);
        assert_eq!(ground_speed_ms(None), None);
    }
}
