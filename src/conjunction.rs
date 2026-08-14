//! Conjunction screening: closest point of approach between track pairs.
//!
//! Under a constant-velocity assumption the separation between two tracks is a
//! quadratic in time, so the closest approach has a closed form and needs no search.
//! With relative position dp and relative velocity dv:
//!
//! ```text
//! t_cpa = -(dp . dv) / (dv . dv)
//! ```
//!
//! clamped to [0, horizon] because a closest approach in the past is not a warning.
//!
//! Horizontal and vertical separation are judged separately, because that is how
//! airspace is actually divided: two aircraft directly above one another are not in
//! conflict if a thousand feet stands between them.

use crate::track::{Track, TrackStore, TrackView};
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct Conjunction {
    pub a: String,
    pub b: String,
    pub label_a: String,
    pub label_b: String,
    /// Seconds from now to closest approach.
    pub t_cpa: f64,
    pub horiz_m: f64,
    pub vert_m: f64,
    /// Combined one-sigma position uncertainty of the pair, metres.
    pub sigma_m: f64,
    /// Closing speed at the moment of detection, m/s.
    pub closing_ms: f64,
}

impl Conjunction {
    /// Severity ordering: soonest first, then tightest.
    fn rank(&self) -> f64 {
        self.t_cpa + self.horiz_m / 1000.0
    }
}

/// Closed-form closest approach between two firm tracks, or `None` if they never come
/// inside the minima within the horizon.
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

    // Cheap reject: even closing head-on at full relative speed, can they reach the
    // minima inside the horizon? This is what keeps an O(n^2) screen affordable.
    let closing = dv.norm();
    if dp.norm() > closing * horizon_s + min_horiz_m {
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
    let horiz = at_cpa.horiz();
    let vert = at_cpa.u.abs();

    if horiz < min_horiz_m && vert < min_vert_m {
        Some(Conjunction {
            a: a.id.clone(),
            b: b.id.clone(),
            label_a: a.label.clone(),
            label_b: b.label.clone(),
            t_cpa: t,
            horiz_m: horiz,
            vert_m: vert,
            // The extrapolated uncertainty, not the uncertainty at the last report: a
            // track that has been coasting for eight seconds deserves a wider ellipse.
            sigma_m: (av.pos_sigma.powi(2) + bv.pos_sigma.powi(2)).sqrt(),
            closing_ms: closing,
        })
    } else {
        None
    }
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
    out.sort_by(|x, y| x.rank().total_cmp(&y.rank()));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geo::{Frame, Geodetic};
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
