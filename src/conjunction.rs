//! Conjunction screening: loss of separation between track pairs.
//!
//! Horizontal and vertical separation are judged separately, because that is how
//! airspace is actually divided: two aircraft directly above one another are not in
//! conflict if a thousand feet stands between them. So the question an alert answers
//! is whether some instant in [0, horizon] has the pair inside *both* minima at once.
//!
//! Under a constant-velocity assumption, with relative position dp and relative
//! velocity dv, each half has a closed form. Vertical separation is linear in time, so
//! it is inside its minimum on one open interval. Horizontally the relative track is a
//! straight line: it passes at a closest distance `d` from the origin, and it is inside
//! the minimum `H` along a chord of half-length `w = √(H² − d²)` centred on that point,
//! which it crosses at the relative ground speed. The pair alerts when the two
//! intervals overlap somewhere in [0, horizon]; where the overlap starts is the time to
//! loss of separation.
//!
//! Whether there is a breach is not decided from the rounded ends of those windows. A
//! window can be narrower than the spacing of doubles at the instant it happens: against
//! a minimum of 1e-15 m, two aircraft passing through each other 1 s from now are inside
//! it for 2e-17 s, both of its rounded ends are 1.0, and the window read as empty. The
//! decision is made instead from the signs of polynomials in the relative state and the
//! minima (see `decide`). Each sign is proven (see `robust`), and a tangent gives exactly
//! zero, which is not a breach. A sign that cannot be proven leaves the pair reported,
//! marked unresolved, rather than dropped.
//!
//! Only the time to loss of separation, once a breach is established, comes from the
//! windows' ends. Those are formed along and across the relative track, with the
//! half-chord `√((H − d)(H + d))`, so that no position, velocity or minimum is squared on
//! its own: a minimum of 1e-300 m squares to 0.
//!
//! The closest approach in three dimensions, `t_cpa` from [`closest_approach_s`], is
//! still reported, but it does not decide anything. Metres of altitude and metres of
//! range are not interchangeable against minima of 305 m and 9260 m, so the instant of
//! least 3-D distance need not lie inside the breach at all.

use crate::geo::Enu;
use crate::track::{Track, TrackStore, TrackView};
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct Conjunction {
    pub a: String,
    pub b: String,
    pub label_a: String,
    pub label_b: String,
    /// Seconds from now until the pair is inside both minima at once; 0 if it already
    /// is. This is what the alert is about. Estimated in floating point once the breach
    /// is established; whether there is one does not depend on it.
    pub t_los: f64,
    /// Whether the breach is proven. `false` means the arithmetic could establish neither
    /// answer, which can happen at extreme scales such as a 1e-300 m minimum, inside the
    /// accepted configuration; it does not mean the pair is safe. The pair is reported
    /// so that it is not silently dropped.
    pub resolved: bool,
    /// Seconds from now to the closest approach in three dimensions. Context only: it
    /// can fall outside the breach.
    pub t_cpa: f64,
    /// Horizontal and vertical separation at `t_cpa`.
    pub horiz_m: f64,
    pub vert_m: f64,
    /// Combined one-sigma position uncertainty of the pair, metres.
    pub sigma_m: f64,
    /// Closing speed at the moment of detection, m/s.
    pub closing_ms: f64,
}

/// Severity ordering: soonest loss of separation first, and tightest only among those
/// that begin together. Adding the two into one score let distance outvote time.
fn by_severity(x: &Conjunction, y: &Conjunction) -> std::cmp::Ordering {
    x.t_los
        .total_cmp(&y.t_los)
        .then(x.horiz_m.total_cmp(&y.horiz_m))
}

/// The breach decision, as signs of polynomials.
mod decide {
    use crate::geo::Enu;
    use crate::robust::{self, Ball, Exact, Inputs, Num, H, PE, PN, PU, T, V, VE, VN, VU};
    use std::cmp::Ordering;

    /// Proven yes, proven no, or not established.
    type Tri = Option<bool>;

    fn and(a: Tri, b: Tri) -> Tri {
        match (a, b) {
            (Some(false), _) | (_, Some(false)) => Some(false),
            (Some(true), Some(true)) => Some(true),
            _ => None,
        }
    }

    fn or(a: Tri, b: Tri) -> Tri {
        match (a, b) {
            (Some(true), _) | (_, Some(true)) => Some(true),
            (Some(false), Some(false)) => Some(false),
            _ => None,
        }
    }

    /// Whether predicate `$f`, reading the inputs `$uses`, has the sign `$want`.
    macro_rules! sign_is {
        ($x:expr, $f:ident, $uses:expr, $want:ident) => {
            robust::sign($x, &$uses, $f::<Ball>, $f::<Exact>).map(|s| s == Ordering::$want)
        };
    }

    // With the vertical velocity made non-negative, the vertical window opens at
    // t1 = −(V + pu)/vu and closes at t2 = (V − pu)/vu. Horizontally, f(t) = |p + vt|² − H²
    // is a convex quadratic with its vertex at t* = −(p·v)/|v|². Each predicate below is
    // one of these quantities times a positive factor that clears its denominators, so
    // it keeps the quantity's sign.

    fn sq<N: Num>(a: &N) -> N {
        a.mul(a)
    }

    fn speed2<N: Num>(x: &[N; 9]) -> N {
        sq(&x[VE]).add(&sq(&x[VN]))
    }

    /// p·v, negative iff the closest approach is still ahead.
    fn approach<N: Num>(x: &[N; 9]) -> N {
        x[PE].mul(&x[VE]).add(&x[PN].mul(&x[VN]))
    }

    /// f(0).
    fn inside_now<N: Num>(x: &[N; 9]) -> N {
        sq(&x[PE]).add(&sq(&x[PN])).sub(&sq(&x[H]))
    }

    /// f(T).
    fn inside_at_horizon<N: Num>(x: &[N; 9]) -> N {
        let e = x[PE].add(&x[VE].mul(&x[T]));
        let n = x[PN].add(&x[VN].mul(&x[T]));
        sq(&e).add(&sq(&n)).sub(&sq(&x[H]))
    }

    /// |v|²·f(t*) = (p × v)² − H²|v|². Negative iff the track passes inside the minimum;
    /// exactly zero for a tangent.
    fn passes_inside<N: Num>(x: &[N; 9]) -> N {
        let cross = x[PE].mul(&x[VN]).sub(&x[PN].mul(&x[VE]));
        sq(&cross).sub(&sq(&x[H]).mul(&speed2(x)))
    }

    /// |v|²·(T − t*).
    fn before_horizon<N: Num>(x: &[N; 9]) -> N {
        approach(x).add(&x[T].mul(&speed2(x)))
    }

    /// vu²·f(t1).
    fn inside_at_opening<N: Num>(x: &[N; 9]) -> N {
        let k = x[V].add(&x[PU]);
        let e = x[PE].mul(&x[VU]).sub(&x[VE].mul(&k));
        let n = x[PN].mul(&x[VU]).sub(&x[VN].mul(&k));
        sq(&e).add(&sq(&n)).sub(&sq(&x[H].mul(&x[VU])))
    }

    /// vu²·f(t2).
    fn inside_at_closing<N: Num>(x: &[N; 9]) -> N {
        let k = x[V].sub(&x[PU]);
        let e = x[PE].mul(&x[VU]).add(&x[VE].mul(&k));
        let n = x[PN].mul(&x[VU]).add(&x[VN].mul(&k));
        sq(&e).add(&sq(&n)).sub(&sq(&x[H].mul(&x[VU])))
    }

    /// |v|²·vu·(t* − t1).
    fn after_opening<N: Num>(x: &[N; 9]) -> N {
        x[V].add(&x[PU])
            .mul(&speed2(x))
            .sub(&approach(x).mul(&x[VU]))
    }

    /// |v|²·vu·(t2 − t*).
    fn before_closing<N: Num>(x: &[N; 9]) -> N {
        x[V].sub(&x[PU])
            .mul(&speed2(x))
            .add(&approach(x).mul(&x[VU]))
    }

    /// vu·(T − t1).
    fn opens_before_horizon<N: Num>(x: &[N; 9]) -> N {
        x[VU].mul(&x[T]).add(&x[V]).add(&x[PU])
    }

    /// Whether the pair is inside both minima at some instant of [0, horizon]: proven
    /// yes, proven no, or `None` if the arithmetic establishes neither.
    ///
    /// The horizontal window, the vertical window and [0, horizon] are intervals, and
    /// intervals on a line share a point iff every two of them do (Helly's theorem in one
    /// dimension). Each of the three pairings asks whether a convex function dips below
    /// zero over an interval. Its minimum there is at the vertex if the vertex lies
    /// inside, and at an end otherwise, so every question is the sign of one polynomial.
    /// The windows are open, because the minima are strict, and a closed end can stand
    /// in for an open one: where an open window's end is strictly inside the other set,
    /// so are the window's points next to it.
    pub(super) fn separation(dp: Enu, dv: Enu, horizon_s: f64, min_h: f64, min_v: f64) -> Tri {
        // The vertical window is unchanged when both its signs flip; flipping makes t1 its
        // opening.
        let (pu, vu) = if dv.u < 0.0 {
            (-dp.u, -dv.u)
        } else {
            (dp.u, dv.u)
        };
        let x: Inputs = [dp.e, dp.n, pu, dv.e, dv.n, vu, min_h, min_v, horizon_s];
        // A minimum that is not positive encloses nothing, and a negative horizon contains
        // no instant. Non-finite state is refused where it enters, and has no breach here.
        if min_h <= 0.0 || min_v <= 0.0 || horizon_s < 0.0 || x.iter().any(|v| !v.is_finite()) {
            return Some(false);
        }

        // The vertical window meets [0, T]: it closes after 0, and opens before T.
        let vertical = if vu == 0.0 {
            Some(pu.abs() < min_v)
        } else {
            and(
                Some(pu < min_v),
                sign_is!(&x, opens_before_horizon, [PU, VU, V, T], Greater),
            )
        };
        if vertical == Some(false) {
            return vertical;
        }
        if dv.e == 0.0 && dv.n == 0.0 {
            // f is constant: horizontally inside throughout, or never.
            return and(vertical, sign_is!(&x, inside_now, [PE, PN, H], Less));
        }
        let crosses = sign_is!(&x, passes_inside, [PE, PN, VE, VN, H], Less);
        if crosses == Some(false) {
            return crosses;
        }
        // The horizontal window meets [0, T]: inside at an end, or at t* between them.
        let within_horizon = or(
            or(
                sign_is!(&x, inside_now, [PE, PN, H], Less),
                sign_is!(&x, inside_at_horizon, [PE, PN, VE, VN, H, T], Less),
            ),
            and(
                and(
                    sign_is!(&x, approach, [PE, PN, VE, VN], Less),
                    sign_is!(&x, before_horizon, [PE, PN, VE, VN, T], Greater),
                ),
                crosses,
            ),
        );
        // The horizontal window meets the vertical one: the same question over [t1, t2].
        let together = if vu == 0.0 {
            crosses
        } else {
            let all = [PE, PN, PU, VE, VN, VU, H, V];
            let order = [PE, PN, PU, VE, VN, VU, V];
            or(
                or(
                    sign_is!(&x, inside_at_opening, all, Less),
                    sign_is!(&x, inside_at_closing, all, Less),
                ),
                and(
                    and(
                        sign_is!(&x, after_opening, order, Greater),
                        sign_is!(&x, before_closing, order, Greater),
                    ),
                    crosses,
                ),
            )
        };
        and(and(vertical, within_horizon), together)
    }
}

/// Where the breach begins, estimated in floating point once `decide::separation` has
/// established it: the later of the two windows' openings, kept inside [0, horizon]. A
/// window narrower than the spacing of doubles collapses here to a point, which is within
/// a rounding of the true instant.
fn entry_estimate(dp: Enu, dv: Enu, horizon_s: f64, min_h: f64, min_v: f64) -> f64 {
    let vertical = if dv.u == 0.0 {
        f64::NEG_INFINITY
    } else {
        ((-min_v - dp.u) / dv.u).min((min_v - dp.u) / dv.u)
    };
    vertical
        .max(horizontal_entry(dp, dv, min_h))
        .clamp(0.0, horizon_s)
}

/// When the relative track enters the horizontal minimum.
fn horizontal_entry(p: Enu, v: Enu, limit: f64) -> f64 {
    let speed = v.e.hypot(v.n);
    if speed == 0.0 || !speed.is_finite() {
        return f64::NEG_INFINITY;
    }
    // Position along the track and distance across it, from unit direction components
    // that cannot exceed 1.
    let (ue, un) = (v.e / speed, v.n / speed);
    let (along, across) = (p.e * ue + p.n * un, (p.e * un - p.n * ue).abs());
    if across >= limit || across.is_nan() {
        // Within a rounding of a tangent: the pair is closest, and inside if at all, here.
        return -along / speed;
    }
    // Half the chord inside the minimum. The product is rounded once before the root, and
    // in binary floating point √(x²) rounds back to exactly x, so a pair approaching along
    // its line of centres enters at the exact instant. When the product is not a normal
    // number the two roots are taken first, which cannot underflow.
    let (gap, reach) = (limit - across, limit + across);
    let chord2 = gap * reach;
    let w = if chord2.is_normal() {
        chord2.sqrt()
    } else {
        gap.sqrt() * reach.sqrt()
    };
    (-w - along) / speed
}

/// Seconds from now to the closest approach in three dimensions, clamped to
/// [0, horizon]: `−(dp·dv)/|dv|²`, and 0 when the relative velocity is exactly zero.
///
/// There is no small-speed threshold. One that treated any squared speed below 1e-9 as
/// parallel reported a pair 1 m apart and closing at 1e-5 m/s as already at its closest,
/// when over a 100,000 s horizon it closes to 0 m. The speed is formed with `hypot` and
/// the division goes through the unit direction, so neither a tiny nor a large relative
/// velocity underflows or overflows on the way; an unbounded quotient clamps to the
/// horizon or to zero.
pub fn closest_approach_s(dp: Enu, dv: Enu, horizon_s: f64) -> f64 {
    let speed = dv.e.hypot(dv.n).hypot(dv.u);
    if speed == 0.0 || !speed.is_finite() {
        return 0.0;
    }
    let along = dp.e * (dv.e / speed) + dp.n * (dv.n / speed) + dp.u * (dv.u / speed);
    (-along / speed).clamp(0.0, horizon_s)
}

/// Loss of separation between two firm tracks, or `None` if they are never inside both
/// minima at once within the horizon.
///
/// Takes views rather than tracks for the geometry, because the two filters are almost
/// never valid for the same instant — one aircraft's last report can be six seconds
/// older than another's. Comparing them where they happen to sit is how phantom
/// conflicts are born; the caller extrapolates both to one instant first.
pub fn pair_cpa(
    a: &Track,
    av: TrackView,
    b: &Track,
    bv: TrackView,
    horizon_s: f64,
    min_horiz_m: f64,
    min_vert_m: f64,
) -> Option<Conjunction> {
    let dp = av.pos - bv.pos;
    let dv = av.vel - bv.vel;

    // Every rejection here is proven. There is deliberately no cheap bound in front: the
    // one that stood here compared a 3-D range to the horizontal minimum, and turned away
    // pairs already inside both minima. A pair the arithmetic cannot decide either way is
    // reported, marked unresolved, rather than dropped.
    let verdict = decide::separation(dp, dv, horizon_s, min_horiz_m, min_vert_m);
    if verdict == Some(false) {
        return None;
    }

    let t = closest_approach_s(dp, dv, horizon_s);
    let at_cpa = dp + dv * t;

    Some(Conjunction {
        a: a.id.clone(),
        b: b.id.clone(),
        label_a: a.label.clone(),
        label_b: b.label.clone(),
        t_los: entry_estimate(dp, dv, horizon_s, min_horiz_m, min_vert_m),
        resolved: verdict == Some(true),
        t_cpa: t,
        horiz_m: at_cpa.horiz(),
        vert_m: at_cpa.u.abs(),
        // The extrapolated uncertainty, not the uncertainty at the last report: a
        // track that has been coasting for eight seconds deserves a wider ellipse.
        sigma_m: (av.pos_sigma.powi(2) + bv.pos_sigma.powi(2)).sqrt(),
        closing_ms: dv.norm(),
    })
}

/// Screen every firm pair in the store, as it stands at `now`. Returns the list ordered
/// by severity.
pub fn screen(
    store: &TrackStore,
    now: Instant,
    horizon_s: f64,
    min_horiz_m: f64,
    min_vert_m: f64,
) -> Vec<Conjunction> {
    // Extrapolate once per track, not once per pair: with 150 tracks that is 150 views
    // instead of 22,350 redundant ones.
    let firm: Vec<(&Track, TrackView)> = store
        .iter()
        .filter(|t| t.is_firm() && !t.on_ground)
        .map(|t| (t, t.view_at(now)))
        .collect();

    let mut out = Vec::new();
    for i in 0..firm.len() {
        for j in (i + 1)..firm.len() {
            let (a, av) = firm[i];
            let (b, bv) = firm[j];
            if let Some(c) = pair_cpa(a, av, b, bv, horizon_s, min_horiz_m, min_vert_m) {
                out.push(c);
            }
        }
    }
    out.sort_by(by_severity);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geo::{Enu, Frame, Geodetic};
    use crate::ingest::Contact;
    use crate::track::TrackStore;
    use std::time::{Duration, Instant};

    fn frame() -> Frame {
        Frame::new(Geodetic {
            lat_deg: 37.9838,
            lon_deg: 23.7275,
            alt_m: 0.0,
        })
    }

    fn plot(id: &str, lat: f64, lon: f64, alt_m: f64) -> Contact {
        Contact {
            id: id.into(),
            label: id.into(),
            kind: "TEST".into(),
            squawk: String::new(),
            emergency: false,
            on_ground: false,
            geo: Geodetic {
                lat_deg: lat,
                lon_deg: lon,
                alt_m,
            },
            gs_kt: None,
            track_deg: None,
            vrate_fpm: None,
            age_s: 0.0,
            source: "TEST",
        }
    }

    /// Two aircraft closing head-on at the same altitude, twelve steps of one second.
    /// `south_reports` is how many of those steps the southern aircraft is heard for,
    /// so a test can leave one track coasting. Returns the store and the instant to
    /// screen it at.
    fn head_on(alt_a: f64, alt_b: f64, south_reports: u64) -> (TrackStore, Instant) {
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let start = Instant::now();
        for step in 0..12 {
            let s = step as f64;
            let mut batch = vec![plot("north", 38.20 - 0.0018 * s, 23.7275, alt_a)];
            if step < south_reports {
                batch.push(plot("south", 37.80 + 0.0018 * s, 23.7275, alt_b));
            }
            store.ingest(&batch, &f, start + Duration::from_secs(step));
        }
        (store, start + Duration::from_secs(11))
    }

    #[test]
    fn detects_a_head_on_conflict() {
        let (store, at) = head_on(10_000.0, 10_000.0, 12);
        let alerts = screen(&store, at, 300.0, 9260.0, 305.0);
        assert_eq!(alerts.len(), 1, "expected exactly one conflicting pair");
        let c = &alerts[0];
        assert!(c.t_cpa > 0.0 && c.t_cpa < 300.0, "t_cpa was {}", c.t_cpa);
        assert!(c.horiz_m < 9260.0);
        assert!(c.closing_ms > 300.0, "closing speed was {}", c.closing_ms);
    }

    #[test]
    fn vertical_separation_clears_the_same_geometry() {
        // Identical horizontal conflict, but 2000 ft apart: legal, and not an alert.
        let (store, at) = head_on(10_000.0, 10_610.0, 12);
        assert!(screen(&store, at, 300.0, 9260.0, 305.0).is_empty());
    }

    #[test]
    fn tracks_heard_at_different_times_are_compared_at_one_instant() {
        // This is what `predict_all` used to guarantee and what views guarantee now.
        // The southern aircraft goes quiet after 6 s, so at screening time the two
        // filters are valid for instants 6 s apart. The screen must still find the
        // conflict, and must widen the pair's uncertainty to admit that half of it is
        // six seconds of extrapolation rather than evidence.
        let (fresh, at_fresh) = head_on(10_000.0, 10_000.0, 12);
        let (coasting, at_coast) = head_on(10_000.0, 10_000.0, 6);

        let a = screen(&fresh, at_fresh, 300.0, 9260.0, 305.0);
        let b = screen(&coasting, at_coast, 300.0, 9260.0, 305.0);

        assert_eq!(b.len(), 1, "a coasting track must still be screened");
        assert!(
            b[0].sigma_m > a[0].sigma_m * 1.5,
            "coasting must widen the pair uncertainty: {} vs {}",
            a[0].sigma_m,
            b[0].sigma_m
        );
    }

    #[test]
    fn parallel_tracks_do_not_alert() {
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let start = Instant::now();
        for step in 0..12 {
            let s = step as f64;
            store.ingest(
                &[
                    plot("left", 38.0 + 0.0018 * s, 23.60, 10_000.0),
                    plot("right", 38.0 + 0.0018 * s, 23.90, 10_000.0),
                ],
                &f,
                start + Duration::from_secs(step),
            );
        }
        assert!(screen(
            &store,
            start + Duration::from_secs(11),
            300.0,
            9260.0,
            305.0
        )
        .is_empty());
    }

    /// Two tracks moving apart, `deg_apart` degrees of latitude between them at t0.
    /// Returns the store and the instant to screen it at.
    fn receding(deg_apart: f64) -> (TrackStore, Instant) {
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let start = Instant::now();
        let half = deg_apart / 2.0;
        for step in 0..12 {
            let s = step as f64;
            store.ingest(
                &[
                    plot("up", 37.9838 + half + 0.0018 * s, 23.7275, 10_000.0),
                    plot("down", 37.9838 - half - 0.0018 * s, 23.7275, 10_000.0),
                ],
                &f,
                start + Duration::from_secs(step),
            );
        }
        (store, start + Duration::from_secs(11))
    }

    #[test]
    fn a_receding_pair_is_not_a_warning() {
        // ~22 km apart and opening. The quadratic minimum is in the past, and clamping
        // to t=0 must not resurrect it as a future conflict.
        let (store, at) = receding(0.20);
        let alerts = screen(&store, at, 300.0, 9260.0, 305.0);
        assert!(
            alerts.is_empty(),
            "separating traffic outside the minima must be quiet"
        );
    }

    #[test]
    fn a_current_violation_alerts_even_while_separating() {
        // Deliberate, not incidental: 5 NM at the same level is a loss of separation
        // now. That it is opening rather than closing does not un-lose it, and an
        // operator still has to see it. t_cpa reads 0 because now is the worst moment.
        let (store, at) = receding(0.01);
        let alerts = screen(&store, at, 300.0, 9260.0, 305.0);
        assert_eq!(alerts.len(), 1);
        assert_eq!(alerts[0].t_cpa, 0.0, "the closest approach is the present");
        assert!(alerts[0].horiz_m < 9260.0);
    }

    fn enu(e: f64, n: f64, u: f64) -> Enu {
        Enu { e, n, u }
    }

    /// Two real tracks, for their identities only.
    fn two_tracks() -> TrackStore {
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        store.ingest(
            &[
                plot("a", 38.0, 23.70, 10_000.0),
                plot("b", 38.0, 23.80, 10_000.0),
            ],
            &frame(),
            Instant::now(),
        );
        store
    }

    fn view((pos, vel): (Enu, Enu)) -> TrackView {
        TrackView {
            pos,
            vel,
            pos_sigma: 0.0,
            coast_s: 0.0,
        }
    }

    /// The views carry the given states verbatim, so a regression pins the numbers it
    /// names rather than whatever a filter happens to converge to.
    fn exact_pair(a: (Enu, Enu), b: (Enu, Enu), horizon_s: f64) -> Option<Conjunction> {
        let store = two_tracks();
        pair_cpa(
            store.get("a").unwrap(),
            view(a),
            store.get("b").unwrap(),
            view(b),
            horizon_s,
            9260.0,
            305.0,
        )
    }

    #[test]
    fn a_breach_away_from_the_3d_closest_approach_alerts() {
        // Audit regression. At t = 60 s the pair is 9000 m apart horizontally and level,
        // inside both minima. The 3-D closest approach falls at ~149 s instead, by which
        // time the vertical gap has opened to ~891 m; judging only that instant missed it.
        // Horizontal is inside on (57.4, 242.6) s, vertical on (29.5, 90.5) s.
        let c = exact_pair(
            (enu(15_000.0, 0.0, -600.0), enu(-100.0, 0.0, 10.0)),
            (enu(0.0, 0.0, 0.0), enu(0.0, 0.0, 0.0)),
            180.0,
        );
        let c = c.expect("simultaneous breach at t = 60 s was not alerted");
        assert!((c.t_los - 57.4).abs() < 1e-9, "t_los was {}", c.t_los);
        assert!(
            c.t_cpa > 90.5,
            "the 3-D CPA ({}) lies outside the breach",
            c.t_cpa
        );
    }

    #[test]
    fn colocation_breaches_any_positive_minimum() {
        // Audit regression. A horizontal minimum of 1e-300 m squared underflows to zero,
        // and "distance² < limit²" became "0 < 0". Exact colocation is inside every
        // positive minimum, stationary or moving together.
        let origin = (enu(0.0, 0.0, 0.0), enu(0.0, 0.0, 0.0));
        let together = (enu(0.0, 0.0, 0.0), enu(120.0, -40.0, 3.0));
        let breach = |a, b, limit: f64| {
            let store = two_tracks();
            let c = pair_cpa(
                store.get("a").unwrap(),
                view(a),
                store.get("b").unwrap(),
                view(b),
                180.0,
                limit,
                305.0,
            );
            let c = c.unwrap_or_else(|| panic!("colocated pair missed at {limit} m"));
            assert!(c.resolved, "colocation at {limit} m was not proven");
            assert_eq!(c.t_los, 0.0);
        };
        // Colocated now and moving apart at 100 m/s: inside for 2·limit/100 s around now,
        // a span no double holds at the smallest limits. The decision does not need it to.
        let through = (enu(0.0, 0.0, 0.0), enu(-100.0, 0.0, 0.0));
        for limit in [
            1e-300,
            1e-200,
            1e-160,
            1e-155,
            f64::MIN_POSITIVE,
            5e-324,
            1e-3,
            9260.0,
        ] {
            breach(origin, origin, limit);
            breach(together, together, limit);
            breach(through, origin, limit);
        }
    }

    /// Pair `a` against a still track at the origin, so the relative state is `a` exactly,
    /// in both orders. The verdict and its estimate must not depend on which comes first.
    /// Returns the estimated onset and whether the breach is proven.
    fn verdict(a: (Enu, Enu), horizon: f64, h: f64, v: f64) -> Option<(f64, bool)> {
        let store = two_tracks();
        let (ta, tb) = (store.get("a").unwrap(), store.get("b").unwrap());
        let still = view((enu(0.0, 0.0, 0.0), enu(0.0, 0.0, 0.0)));
        let one = pair_cpa(ta, view(a), tb, still, horizon, h, v).map(|c| (c.t_los, c.resolved));
        let other = pair_cpa(tb, still, ta, view(a), horizon, h, v).map(|c| (c.t_los, c.resolved));
        assert_eq!(
            one, other,
            "reversing the pair changed the verdict for {a:?}"
        );
        one
    }

    #[test]
    fn a_breach_narrower_than_the_spacing_of_doubles_is_found() {
        // Audit regression, the exact witness. 100 m apart, closing at 100 m/s, level,
        // against a 1e-15 m minimum: exactly colocated at t = 1 s, inside for 2e-17 s.
        // Both rounded ends of that window are 1.0, and reading the window from them
        // found it empty.
        assert_eq!((100.0 - 1e-15) / 100.0, (100.0 + 1e-15) / 100.0);
        let head_on = |pn: f64| (enu(-100.0, pn, 0.0), enu(100.0, 0.0, 0.0));
        let (t_los, proven) = verdict(head_on(0.0), 180.0, 1e-15, 305.0).expect("missed");
        assert!(proven && (t_los - 1.0).abs() < 1e-15, "{t_los} {proven}");

        // At the same scale, a tangent is not a breach, one ulp outside is not, and one
        // ulp inside is.
        let edge = 1e-15f64;
        assert_eq!(verdict(head_on(edge), 180.0, edge, 305.0), None, "tangent");
        assert_eq!(
            verdict(head_on(edge.next_up()), 180.0, edge, 305.0),
            None,
            "miss"
        );
        assert!(verdict(head_on(edge.next_down()), 180.0, edge, 305.0).is_some_and(|v| v.1));

        // The horizon boundary at the same scale: colocated exactly at a 1 s horizon is a
        // breach, and a horizon one ulp shorter ends before the window opens.
        assert!(verdict(head_on(0.0), 1.0, 1e-15, 305.0).is_some_and(|v| v.1));
        assert_eq!(verdict(head_on(0.0), 1f64.next_down(), 1e-15, 305.0), None);
    }

    #[test]
    fn a_tangent_is_not_a_breach_in_any_direction() {
        // At an ordinary scale, in a direction whose unit vector no double holds exactly:
        // with v = (3, 4), the track from (740, −680) passes exactly 1000 m from the
        // origin, 20 s from now. The across-track distance, computed through 0.6 and 0.8,
        // lands either side of 1000 by rounding; the decision may not.
        let pass = (enu(740.0, -680.0, 0.0), enu(3.0, 4.0, 0.0));
        assert_eq!(verdict(pass, 180.0, 1000.0, 305.0), None);
        let (t_los, proven) = verdict(pass, 180.0, 1000f64.next_up(), 305.0).expect("missed");
        assert!(proven && (t_los - 20.0).abs() < 1e-3, "{t_los}");
        // An ordinary crossing, and its onset: 50 km apart, closing at 500 m/s.
        let crossing = (enu(-30_000.0, -40_000.0, 0.0), enu(300.0, 400.0, 0.0));
        let (t_los, proven) = verdict(crossing, 180.0, 9260.0, 305.0).expect("missed");
        assert!(proven && (t_los - 81.48).abs() < 1e-9, "{t_los}");
    }

    #[test]
    fn narrow_windows_are_compared_without_rounding_them() {
        let narrow = |pu: f64, vu: f64| (enu(-100.0, 0.0, pu), enu(100.0, 0.0, vu));
        // The 2e-17 s horizontal window inside a vertical window from 0.5 s to 1.5 s, and
        // before one from 1.5 s to 2.5 s.
        let (t_los, proven) = verdict(narrow(-1000.0, 1000.0), 180.0, 1e-15, 500.0).unwrap();
        assert!(proven && (t_los - 1.0).abs() < 1e-15, "{t_los}");
        assert_eq!(verdict(narrow(-2000.0, 1000.0), 180.0, 1e-15, 500.0), None);

        // Both windows narrower than the spacing of doubles at 1 s, and their centres
        // closer together than that spacing: horizontally 1 s ± H/3, vertically
        // (1 + 2^-51/3) s ± V/3. They overlap iff H + V > 2^-51, about 4.4e-16, and both
        // round to the single instant 1.0 either way.
        let (pu, vu) = (-3f64.next_up(), 3.0);
        let both = (enu(-3.0, 0.0, pu), enu(3.0, 0.0, vu));
        assert!(
            verdict(both, 180.0, 3e-16, 3e-16).is_some_and(|v| v.1),
            "overlap"
        );
        assert_eq!(verdict(both, 180.0, 2e-16, 2e-16), None, "gap");
    }

    #[test]
    fn an_undecidable_pair_is_reported_not_dropped() {
        // Far outside the envelope, through the public function: 1e300 m apart, closing at
        // 1e300 m/s. The exact products leave the range where they are exact, so neither
        // answer can be proven, and the pair is reported, marked.
        let far = (enu(-1e300, 0.0, 0.0), enu(1e300, 0.0, 0.0));
        let (t_los, proven) = verdict(far, 180.0, 9260.0, 305.0).expect("dropped");
        assert!(!proven && t_los.is_finite(), "{t_los} {proven}");
    }

    #[test]
    fn an_approach_along_the_line_of_centres_enters_at_the_exact_instant() {
        // Separation re-audit controls, kept here after a change of formula broke them:
        // √9260 · √9260 is not 9260 in floating point, and the entry time came out as
        // 1.00000000000002 s. 100 m outside the minimum, closing at 100 m/s, the pair
        // enters at 1 s exactly; on the boundary and closing, it enters now.
        let still = (enu(0.0, 0.0, 0.0), enu(0.0, 0.0, 0.0));
        let inbound = |e: f64, u: f64, vu: f64| (enu(e, 0.0, u), enu(-100.0, 0.0, vu));
        let t_los = |a, horizon| exact_pair(a, still, horizon).map(|c| c.t_los);
        assert_eq!(t_los(inbound(9_360.0, 0.0, 0.0), 1.001), Some(1.0));
        assert_eq!(t_los(inbound(9_360.0, 0.0, 0.0), 1.0), None);
        assert_eq!(t_los(inbound(9_260.0, 0.0, 0.0), 180.0), Some(0.0));
        // Vertically inside only until 1.01 s: a 10 ms overlap, starting at exactly 1 s.
        assert_eq!(t_los(inbound(9_360.0, 303.99, 1.0), 180.0), Some(1.0));
    }

    #[test]
    fn a_pair_already_inside_both_minima_is_not_rejected_cheaply() {
        // Audit regression. 9259 m horizontal and 300 m vertical at equal velocity: inside
        // both minima now and for the whole horizon. The 3-D range, ~9263.9 m, exceeds
        // the horizontal minimum, and a 3-D range is what the old reject compared to it.
        let c = exact_pair(
            (enu(9_259.0, 0.0, 300.0), enu(230.0, 0.0, 0.0)),
            (enu(0.0, 0.0, 0.0), enu(230.0, 0.0, 0.0)),
            180.0,
        );
        let c = c.expect("a pair inside both minima now was not alerted");
        assert_eq!(c.t_los, 0.0);
    }

    fn alert(t_los: f64, horiz_m: f64) -> Conjunction {
        Conjunction {
            a: String::new(),
            b: String::new(),
            label_a: String::new(),
            label_b: String::new(),
            t_los,
            resolved: true,
            t_cpa: t_los,
            horiz_m,
            vert_m: 0.0,
            sigma_m: 0.0,
            closing_ms: 0.0,
        }
    }

    #[test]
    fn an_earlier_loss_of_separation_always_ranks_first() {
        // Audit regression. Onset plus CPA distance / 1000 scored a breach already under
        // way at 9259 m as 9.259, and one starting in 1 s at 0 m as 1.0, so the later
        // breach was listed first.
        let mut alerts = [alert(1.0, 0.0), alert(0.0, 9_259.0)];
        alerts.sort_by(by_severity);
        assert_eq!(alerts[0].t_los, 0.0, "the breach under way must lead");

        // Distance decides only between breaches that begin together.
        let mut tied = [alert(5.0, 800.0), alert(5.0, 200.0)];
        tied.sort_by(by_severity);
        assert_eq!(tied[0].horiz_m, 200.0);
    }

    /// A deterministic generator, so a failing case can be named and replayed.
    struct Lcg(u64);

    impl Lcg {
        fn range(&mut self, lo: f64, hi: f64) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            lo + (hi - lo) * ((self.0 >> 11) as f64 / (1u64 << 53) as f64)
        }
    }

    #[test]
    fn the_screen_agrees_with_a_sampled_oracle() {
        // Properties over generated geometry rather than cases someone thought of:
        //   * any sampled instant inside both minima is alerted;
        //   * an alert's onset is real: on or inside both minima at t_los, and no sampled
        //     instant before it is inside both.
        // The oracle only evaluates positions and shares no algebra with the screen. A
        // breach shorter than the sampling step is invisible to it, so this cannot rule
        // out short misses; it can catch the screen disagreeing with plain geometry.
        const H: f64 = 9260.0;
        const V: f64 = 305.0;
        const STEP: f64 = 0.05;
        let store = two_tracks();
        let (ta, tb) = (store.get("a").unwrap(), store.get("b").unwrap());
        let still = (Enu::default(), Enu::default());
        let mut rng = Lcg(0x5eed);
        let (mut alerted, mut quiet) = (0, 0);

        for case in 0..4000 {
            // Every choice is drawn, not keyed to the case number, so each kind of motion
            // meets each kind of placement and both horizons.
            let horizon = if rng.range(0.0, 1.0) < 0.05 {
                0.0
            } else {
                180.0
            };
            let mut dv = enu(
                rng.range(-300.0, 300.0),
                rng.range(-300.0, 300.0),
                rng.range(-40.0, 40.0),
            );
            match rng.range(0.0, 8.0) as u32 {
                0 => dv.u = 0.0,
                1 => (dv.e, dv.n) = (0.0, 0.0),
                // Equal velocities: formation flight, or two airliners at one speed.
                2 => dv = Enu::default(),
                // Barely closing, under 5 cm/s.
                3 => dv = dv * 1e-4,
                _ => {}
            }
            // Half the cases are aimed: a point near the cylinder at some instant, walked
            // back to t = 0, so windows are short and the 3-D closest approach usually
            // lies somewhere else. Half of those sit within 30 m of the rim, where wall
            // meets lid; that corner is where a 3-D shortcut and the true test disagree.
            let dp = if rng.range(0.0, 1.0) < 0.5 {
                let tc = rng.range(0.0, 200.0);
                let th = rng.range(0.0, std::f64::consts::TAU);
                let (r, u) = if rng.range(0.0, 1.0) < 0.5 {
                    let side = if rng.range(0.0, 1.0) < 0.5 { -1.0 } else { 1.0 };
                    (
                        H + rng.range(-30.0, 30.0),
                        side * (V + rng.range(-30.0, 30.0)),
                    )
                } else {
                    (rng.range(0.0, 1.2 * H), rng.range(-1.3 * V, 1.3 * V))
                };
                enu(r * th.cos(), r * th.sin(), u) - dv * tc
            } else {
                enu(
                    rng.range(-30e3, 30e3),
                    rng.range(-30e3, 30e3),
                    rng.range(-1500.0, 1500.0),
                )
            };

            let at = |t: f64| dp + dv * t;
            let inside = |p: Enu| p.horiz() < H && p.u.abs() < V;
            let samples = (horizon / STEP) as usize;
            let first = (0..=samples)
                .map(|k| k as f64 * STEP)
                .find(|&t| inside(at(t)));

            match (
                first,
                pair_cpa(ta, view((dp, dv)), tb, view(still), horizon, H, V),
            ) {
                (Some(t), None) => {
                    panic!("case {case}: inside both at {t} s, no alert. dp {dp:?} dv {dv:?}")
                }
                (first, Some(c)) => {
                    let p = at(c.t_los);
                    assert!(
                        (0.0..=horizon).contains(&c.t_los),
                        "case {case}: onset {} outside the horizon",
                        c.t_los
                    );
                    assert!(
                        p.horiz() <= H + 1e-6 && p.u.abs() <= V + 1e-6,
                        "case {case}: not in conflict at its own onset {} s: {p:?}",
                        c.t_los
                    );
                    if let Some(t) = first {
                        assert!(
                            t >= c.t_los - 1e-6,
                            "case {case}: inside both at {t} s, before onset {} s",
                            c.t_los
                        );
                    }
                    alerted += 1;
                }
                (None, None) => quiet += 1,
            }
        }
        // A generator that stopped producing one side would pass vacuously.
        println!("alerted {alerted}, quiet {quiet}");
        assert!(
            alerted > 500 && quiet > 500,
            "alerted {alerted}, quiet {quiet}"
        );
    }

    #[test]
    fn unfirm_tracks_are_excluded() {
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let now = Instant::now();
        store.ingest(
            &[
                plot("a", 38.0, 23.7275, 10_000.0),
                plot("b", 38.001, 23.7275, 10_000.0),
            ],
            &f,
            now,
        );
        assert!(screen(&store, now, 300.0, 9260.0, 305.0).is_empty());
    }
}
