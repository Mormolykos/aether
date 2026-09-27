//! Conjunction screening: loss of separation between track pairs.
//!
//! Horizontal and vertical separation are judged separately, because that is how
//! airspace is actually divided: two aircraft directly above one another are not in
//! conflict if a thousand feet stands between them. So the question an alert answers
//! is whether some instant in [0, horizon] has the pair inside *both* minima at once.
//!
//! Under a constant-velocity assumption, with relative position dp and relative
//! velocity dv, each half has a closed form. Vertical separation is linear in time, so
//! it is inside its minimum on one open interval. Horizontal separation squared is a
//! convex quadratic,
//!
//! ```text
//! |dv_h|² t² + 2 (dp_h . dv_h) t + |dp_h|² - H² < 0
//! ```
//!
//! inside its minimum on the open interval between the roots. The pair alerts when the
//! two intervals overlap somewhere in [0, horizon]; where the overlap starts is the time
//! to loss of separation.
//!
//! The closest approach in three dimensions, `t_cpa = -(dp . dv) / (dv . dv)` clamped
//! to [0, horizon], is still reported, but it does not decide anything. Metres of
//! altitude and metres of range are not interchangeable against minima of 305 m and
//! 9260 m, so the instant of least 3-D distance need not lie inside the breach at all.

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
    /// is. This is what the alert is about.
    pub t_los: f64,
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

/// The open interval between two roots, or `None` if it is empty. A NaN end gives no
/// interval, as NaN gave no alert before this function existed; what non-finite state
/// should raise instead is a question for where it enters, not for the geometry.
fn between(t1: f64, t2: f64) -> Option<(f64, f64)> {
    let (lo, hi) = if t1 <= t2 { (t1, t2) } else { (t2, t1) };
    (lo < hi).then_some((lo, hi))
}

/// When `|p + v t| < limit` along one axis: linear in `t`, so one open interval.
fn inside_vertical(p: f64, v: f64, limit: f64) -> Option<(f64, f64)> {
    if v == 0.0 {
        return (p.abs() < limit).then_some((f64::NEG_INFINITY, f64::INFINITY));
    }
    between((-limit - p) / v, (limit - p) / v)
}

/// When `|p_h + v_h t| < limit` in the horizontal plane: `a t² + b t + c < 0` with
/// `a >= 0`, so one open interval between the roots.
fn inside_horizontal(p: Enu, v: Enu, limit: f64) -> Option<(f64, f64)> {
    let a = v.e * v.e + v.n * v.n;
    let b = 2.0 * (p.e * v.e + p.n * v.n);
    let c = p.e * p.e + p.n * p.n - limit * limit;
    if a == 0.0 {
        return (c < 0.0).then_some((f64::NEG_INFINITY, f64::INFINITY));
    }
    let disc = b * b - 4.0 * a * c;
    // A tangent only touches the minimum; it never goes inside it.
    if disc.is_nan() || disc <= 0.0 {
        return None;
    }
    // The stable form. `-b ± sqrt(disc)` loses the near root to cancellation when b²
    // dwarfs 4ac, which is a fast pair at long range; q never has that subtraction.
    let q = -0.5 * (b + b.signum() * disc.sqrt());
    between(q / a, c / q)
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

    // Every rejection here comes from the interval algebra, not from a bound. There is
    // deliberately no cheap bound in front: the one that stood here compared a 3-D range
    // to the horizontal minimum, and turned away pairs already inside both minima. The
    // algebra is exact under constant velocity; its evaluation is floating point, so a
    // pair within rounding of a limit can land either side of it, and inputs far outside
    // airspace scales can overflow the squared terms below.
    let (v_lo, v_hi) = inside_vertical(dp.u, dv.u, min_vert_m)?;
    let (h_lo, h_hi) = inside_horizontal(dp, dv, min_horiz_m)?;
    // Inside both at once on the open interval (lo, hi), which has to meet the closed
    // interval [0, horizon]. An empty interval stays empty under the intersection.
    let (lo, hi) = (v_lo.max(h_lo), v_hi.min(h_hi));
    let meets_horizon = lo < hi && lo < horizon_s && hi > 0.0;
    if !meets_horizon {
        return None;
    }

    let dvv = dv.e * dv.e + dv.n * dv.n + dv.u * dv.u;
    // Parallel tracks: separation never changes, so evaluate it now.
    let t = if dvv < 1e-9 {
        0.0
    } else {
        let t = -(dp.e * dv.e + dp.n * dv.n + dp.u * dv.u) / dvv;
        t.clamp(0.0, horizon_s)
    };
    let at_cpa = dp + dv * t;

    Some(Conjunction {
        a: a.id.clone(),
        b: b.id.clone(),
        label_a: a.label.clone(),
        label_b: b.label.clone(),
        t_los: lo.max(0.0),
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
