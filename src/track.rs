//! Track store: association, state estimation, ageing.
//!
//! A contact is a plot; a track is a hypothesis about an object that persists across
//! plots. ADS-B hands over a unique ICAO address, so association is by identity rather
//! than nearest-neighbour gating. The interesting work is what happens after: placing
//! each measurement at the moment it describes, rejecting reports that cannot be true,
//! and dropping tracks whose evidence has gone stale.
//!
//! # Temporal invariant
//!
//! A filter's state is valid for one instant and one only. Here that instant is the
//! observation `valid_age` seconds before the batch that arrived at `valid_rx`, and it
//! is deliberately **never materialised** — only differences of it are ever computed.
//! Both terms are things that really happened: an arrival the process witnessed, and an
//! age the sensor reported. Nothing subtracts a duration from an `Instant`, so nothing
//! can underflow, and no clock the process does not own is ever trusted.
//!
//! The validity time advances in exactly one place: `Track::update`, when a measurement
//! describes a moment later than the one the filter already holds. The operator picture
//! and the conjunction screen do **not** advance it — they take a `TrackView`, which
//! extrapolates a copy. That separation is what makes it safe to place a measurement in
//! its own past relative to the display: the display never wrote to the filter.

use crate::geo::{Enu, Frame, Geodetic};
use crate::ingest::Contact;
use crate::kalman::Kf1D;
use std::collections::HashMap;
use std::time::{Duration, Instant};

#[derive(Clone, Debug)]
pub struct Track {
    pub id: String,
    pub label: String,
    pub kind: String,
    pub squawk: String,
    pub emergency: bool,
    pub on_ground: bool,
    pub source: &'static str,
    pub geo: Geodetic,

    /// East / North / Up estimators.
    filt: [Kf1D; 3],

    /// The filter is valid for the observation made `valid_age` seconds before the
    /// batch that arrived at `valid_rx`. See the module-level temporal invariant.
    valid_rx: Instant,
    valid_age: f64,

    /// The same pair for the last *accepted* measurement. A plot that coasts the filter
    /// forward but then fails the gate moves `valid_*` and leaves these alone, so
    /// ageing and the staleness column keep meaning "evidence", not "traffic".
    accept_rx: Instant,
    accept_age: f64,

    pub first_seen: Instant,
    pub updates: u64,
    /// Plots refused by the innovation gate.
    pub rejected: u64,
    /// Plots describing a moment the filter had already passed: re-served snapshots and
    /// out-of-order reports.
    pub superseded: u64,
    /// Normalised innovation squared of the last accepted update.
    pub last_nis: f64,
}

/// A track's estimate extrapolated to some instant, without disturbing the track.
///
/// Everything that wants to know where an aircraft is *now* takes one of these. The
/// filter itself stays parked at the moment it was last given evidence for.
#[derive(Clone, Copy, Debug)]
pub struct TrackView {
    pub pos: Enu,
    pub vel: Enu,
    pub pos_sigma: f64,
    /// Seconds of coast between the filter's validity time and the viewed instant.
    pub coast_s: f64,
}

/// What became of one plot.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Outcome {
    Initiated,
    Updated,
    /// Refused by the innovation gate.
    Gated,
    /// Describes a moment the filter had already passed.
    Superseded,
}

impl Track {
    fn new(c: &Contact, pos: Enu, received: Instant, q: f64, r: f64) -> Self {
        // Seed velocity from the reported heading and ground speed when present. It is
        // only a prior: the covariance is loose enough that two updates overrule it.
        let (ve, vn) = match (c.gs_kt, c.track_deg) {
            (Some(gs), Some(trk)) => {
                let speed = gs * crate::geo::KT_TO_MS;
                let rad = trk.to_radians();
                (speed * rad.sin(), speed * rad.cos())
            }
            _ => (0.0, 0.0),
        };
        let vu = c.vrate_fpm.unwrap_or(0.0) * crate::geo::FPM_TO_MS;

        Self {
            id: c.id.clone(),
            label: c.label.clone(),
            kind: c.kind.clone(),
            squawk: c.squawk.clone(),
            emergency: c.emergency,
            on_ground: c.on_ground,
            source: c.source,
            geo: c.geo,
            filt: [
                Kf1D::new(pos.e, ve, q, r),
                Kf1D::new(pos.n, vn, q, r),
                Kf1D::new(pos.u, vu, q, r),
            ],
            // The track is initiated at the moment its first plot describes, expressed
            // as the arrival plus the age the sensor reported. No subtraction, so a
            // process started seconds after boot has nothing to underflow.
            valid_rx: received,
            valid_age: c.age_s,
            accept_rx: received,
            accept_age: c.age_s,
            first_seen: received,
            updates: 1,
            rejected: 0,
            superseded: 0,
            last_nis: 0.0,
        }
    }

    /// Seconds of coast between the observation the filter is valid for and `t`.
    ///
    /// This is the single expression the whole temporal design rests on. It is the sum
    /// of an elapsed duration between two witnessed instants and a reported age, so it
    /// is non-negative for any `t` at or after the batch that last touched this track,
    /// and it never reconstructs an absolute observation time.
    pub fn filter_age_at(&self, t: Instant) -> f64 {
        t.saturating_duration_since(self.valid_rx).as_secs_f64() + self.valid_age
    }

    /// The estimate extrapolated to `t`, leaving the filter untouched.
    pub fn view_at(&self, t: Instant) -> TrackView {
        let coast_s = self.filter_age_at(t);
        // `Kf1D` is `Copy`, so this is a register-width snapshot, not an allocation.
        let mut f = self.filt;
        if coast_s > 0.0 {
            for k in &mut f {
                k.predict(coast_s);
            }
        }
        TrackView {
            pos: Enu {
                e: f[0].x,
                n: f[1].x,
                u: f[2].x,
            },
            vel: Enu {
                e: f[0].v,
                n: f[1].v,
                u: f[2].v,
            },
            pos_sigma: f.iter().map(|k| k.pos_sigma()).fold(0.0, f64::max),
            coast_s,
        }
    }

    /// The estimate at the moment the filter is actually valid for, unextrapolated.
    pub fn position(&self) -> Enu {
        Enu {
            e: self.filt[0].x,
            n: self.filt[1].x,
            u: self.filt[2].x,
        }
    }

    pub fn velocity(&self) -> Enu {
        Enu {
            e: self.filt[0].v,
            n: self.filt[1].v,
            u: self.filt[2].v,
        }
    }

    /// Ground speed from the estimator, knots.
    pub fn speed_kt(&self) -> f64 {
        self.velocity().horiz() * crate::geo::MS_TO_KT
    }

    /// Estimated heading over the ground, degrees true.
    pub fn heading_deg(&self) -> f64 {
        let v = self.velocity();
        let h = v.e.atan2(v.n).to_degrees();
        if h < 0.0 {
            h + 360.0
        } else {
            h
        }
    }

    /// Vertical rate, feet per minute.
    pub fn vrate_fpm(&self) -> f64 {
        self.filt[2].v / crate::geo::FPM_TO_MS
    }

    /// Worst per-axis one-sigma position uncertainty, metres.
    pub fn pos_sigma(&self) -> f64 {
        self.filt.iter().map(|f| f.pos_sigma()).fold(0.0, f64::max)
    }

    /// Seconds since the last accepted measurement was *observed* — not since it was
    /// received. A feed that hands us a six-second-old position has given us a
    /// six-second-old track, and the ageing logic should say so.
    pub fn staleness(&self, now: Instant) -> f64 {
        now.saturating_duration_since(self.accept_rx).as_secs_f64() + self.accept_age
    }

    /// A track is firm once enough reports have agreed with it. Only firm tracks are
    /// screened for conjunctions: a two-plot track has a velocity that is mostly prior.
    pub fn is_firm(&self) -> bool {
        self.updates >= 4
    }

    fn update(&mut self, c: &Contact, pos: Enu, received: Instant, gate_sigma: f64) -> Outcome {
        // Elapsed time between the observation the filter holds and the observation
        // this plot describes. Both sides are relative, so this is ordinary f64
        // arithmetic on two witnessed instants and two reported ages.
        let dt = self.filter_age_at(received) - c.age_s;

        if dt <= 0.0 {
            // The plot describes a moment at or before the one the filter already
            // holds. Because the validity time is advanced only by measurements, this
            // is a real property of the data — a re-served snapshot or an out-of-order
            // report — and not an artefact of the display clock. A re-served snapshot
            // lands at dt == 0 exactly, since its arrival and its reported age advance
            // together, so deduplication falls out of the same comparison.
            //
            // v1 drops these. Folding a late measurement back into a filter that has
            // moved past it is retrodiction, and pretending to do it by applying it at
            // the wrong time would be worse than admitting we do not.
            self.superseded += 1;
            return Outcome::Superseded;
        }

        // Coasting to the measurement's own moment is correct whether or not the
        // measurement survives the gate, so the validity time moves first and stays
        // moved. A gated plot leaves a coasted track, which is what it should leave.
        for f in &mut self.filt {
            f.predict(dt);
        }
        self.valid_rx = received;
        self.valid_age = c.age_s;

        // Gate on all three axes before touching any of them: a partially applied
        // update on a bad plot is worse than a rejected one.
        let innovs = [
            self.filt[0].innovation(pos.e),
            self.filt[1].innovation(pos.n),
            self.filt[2].innovation(pos.u),
        ];
        if self.is_firm() && innovs.iter().any(|i| i.sigma() > gate_sigma) {
            self.rejected += 1;
            return Outcome::Gated;
        }

        self.filt[0].update(pos.e);
        self.filt[1].update(pos.n);
        self.filt[2].update(pos.u);

        self.last_nis = innovs.iter().map(|i| i.nis()).sum::<f64>() / 3.0;
        self.updates += 1;
        self.accept_rx = received;
        self.accept_age = c.age_s;
        self.geo = c.geo;
        self.label = c.label.clone();
        self.squawk = c.squawk.clone();
        self.emergency = c.emergency;
        self.on_ground = c.on_ground;
        if !c.kind.is_empty() {
            self.kind = c.kind.clone();
        }
        Outcome::Updated
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub struct IngestReport {
    pub initiated: usize,
    pub updated: usize,
    pub gated: usize,
    pub superseded: usize,
}

pub struct TrackStore {
    tracks: HashMap<String, Track>,
    q: f64,
    r: f64,
    gate_sigma: f64,
    pub total_initiated: u64,
    pub total_dropped: u64,
    pub total_gated: u64,
    pub total_superseded: u64,
}

impl TrackStore {
    pub fn new(process_noise: f64, meas_var: f64, gate_sigma: f64) -> Self {
        Self {
            tracks: HashMap::new(),
            q: process_noise,
            r: meas_var,
            gate_sigma,
            total_initiated: 0,
            total_dropped: 0,
            total_gated: 0,
            total_superseded: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.tracks.len()
    }

    pub fn is_empty(&self) -> bool {
        self.tracks.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &Track> {
        self.tracks.values()
    }

    pub fn get(&self, id: &str) -> Option<&Track> {
        self.tracks.get(id)
    }

    /// Fold one batch in. `received` is when the batch landed; each contact carries how
    /// old it already was at that moment.
    pub fn ingest(
        &mut self,
        contacts: &[Contact],
        frame: &Frame,
        received: Instant,
    ) -> IngestReport {
        let mut report = IngestReport::default();
        for c in contacts {
            let pos = frame.to_enu(c.geo);
            let outcome = match self.tracks.get_mut(&c.id) {
                Some(track) => track.update(c, pos, received, self.gate_sigma),
                None => {
                    self.tracks
                        .insert(c.id.clone(), Track::new(c, pos, received, self.q, self.r));
                    Outcome::Initiated
                }
            };
            match outcome {
                Outcome::Initiated => {
                    report.initiated += 1;
                    self.total_initiated += 1;
                }
                Outcome::Updated => report.updated += 1,
                Outcome::Gated => {
                    report.gated += 1;
                    self.total_gated += 1;
                }
                Outcome::Superseded => {
                    report.superseded += 1;
                    self.total_superseded += 1;
                }
            }
        }
        report
    }

    /// Drop tracks whose last accepted *observation* is older than `timeout`. Returns
    /// how many went.
    pub fn prune(&mut self, now: Instant, timeout: Duration) -> usize {
        let before = self.tracks.len();
        let limit = timeout.as_secs_f64();
        self.tracks.retain(|_, t| t.staleness(now) < limit);
        let dropped = before - self.tracks.len();
        self.total_dropped += dropped as u64;
        dropped
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frame() -> Frame {
        Frame::new(Geodetic {
            lat_deg: 37.9838,
            lon_deg: 23.7275,
            alt_m: 0.0,
        })
    }

    fn contact(id: &str, lat: f64, lon: f64, alt_m: f64) -> Contact {
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
            gs_kt: Some(400.0),
            track_deg: Some(0.0),
            vrate_fpm: Some(0.0),
            age_s: 0.0,
            source: "TEST",
        }
    }

    #[test]
    fn same_id_updates_rather_than_duplicates() {
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let now = Instant::now();
        store.ingest(&[contact("aaa", 38.0, 23.7, 10_000.0)], &f, now);
        store.ingest(
            &[contact("aaa", 38.01, 23.7, 10_000.0)],
            &f,
            now + Duration::from_secs(1),
        );
        assert_eq!(store.len(), 1);
        assert_eq!(store.get("aaa").unwrap().updates, 2);
    }

    #[test]
    fn tracks_a_northbound_target_and_recovers_its_velocity() {
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let start = Instant::now();
        // ~0.0018 deg of latitude per second is about 200 m/s northbound.
        for step in 0..25 {
            let lat = 38.0 + 0.0018 * step as f64;
            store.ingest(
                &[contact("bbb", lat, 23.7, 10_000.0)],
                &f,
                start + Duration::from_secs(step),
            );
        }
        let t = store.get("bbb").unwrap();
        assert!(t.is_firm());
        assert!(
            (t.velocity().n - 200.0).abs() < 25.0,
            "vn was {}",
            t.velocity().n
        );
        assert!(
            t.heading_deg() < 5.0 || t.heading_deg() > 355.0,
            "hdg {}",
            t.heading_deg()
        );
    }

    #[test]
    fn a_teleporting_plot_is_gated_out() {
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let start = Instant::now();
        for step in 0..10 {
            store.ingest(
                &[contact("ccc", 38.0, 23.7, 10_000.0)],
                &f,
                start + Duration::from_secs(step),
            );
        }
        let held = store.get("ccc").unwrap().position();
        // Same aircraft, suddenly 600 km away one second later. Physically impossible.
        let report = store.ingest(
            &[contact("ccc", 43.4, 23.7, 10_000.0)],
            &f,
            start + Duration::from_secs(11),
        );
        assert_eq!(report.gated, 1);
        let after = store.get("ccc").unwrap();
        assert_eq!(after.rejected, 1);
        assert!(
            (after.position() - held).norm() < 5_000.0,
            "the jump was absorbed"
        );
    }

    #[test]
    fn stale_tracks_are_dropped() {
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let now = Instant::now();
        store.ingest(&[contact("ddd", 38.0, 23.7, 10_000.0)], &f, now);
        assert_eq!(
            store.prune(now + Duration::from_secs(10), Duration::from_secs(45)),
            0
        );
        assert_eq!(
            store.prune(now + Duration::from_secs(60), Duration::from_secs(45)),
            1
        );
        assert!(store.is_empty());
        assert_eq!(store.total_dropped, 1);
    }

    #[test]
    fn coasting_widens_uncertainty_but_keeps_the_track() {
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let start = Instant::now();
        for step in 0..10 {
            store.ingest(
                &[contact("eee", 38.0 + 0.0018 * step as f64, 23.7, 10_000.0)],
                &f,
                start + Duration::from_secs(step),
            );
        }
        let t = store.get("eee").unwrap();
        let tight = t.pos_sigma();
        let loose = t.view_at(start + Duration::from_secs(40)).pos_sigma;
        assert!(loose > tight * 2.0, "sigma {tight} -> {loose}");
    }

    // --- temporal invariant ---------------------------------------------------------

    /// Same contact, but the feed says the position was already `age` seconds old.
    fn aged(id: &str, lat: f64, lon: f64, alt_m: f64, age: f64) -> Contact {
        Contact {
            age_s: age,
            ..contact(id, lat, lon, alt_m)
        }
    }

    #[test]
    fn viewing_a_track_never_advances_its_filter() {
        // The defect the whole redesign exists to prevent: the picture must be a
        // reader. If drawing the screen moves the filter, the next measurement is
        // compared against a state from the future and gated for being honest.
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let start = Instant::now();
        store.ingest(&[contact("aaa", 38.0, 23.7, 10_000.0)], &f, start);

        let before = store.get("aaa").unwrap().position();
        for step in 1..=30 {
            let _ = store
                .get("aaa")
                .unwrap()
                .view_at(start + Duration::from_secs(step));
        }
        let after = store.get("aaa").unwrap();
        assert_eq!(after.position(), before, "a view mutated the filter");
        assert_eq!(after.filter_age_at(start), 0.0, "validity time moved");
    }

    /// Reports landing on a fixed poll cadence, each carrying its own age, with that
    /// age jittering from poll to poll. Yields `(arrival_s, age_s, observed_at_s)`.
    ///
    /// The shape is taken from the live Athens feed, measured over two consecutive
    /// polls: median reported age 0.31 s, p90 3.97 s, and a per-aircraft change in age
    /// between polls spanning −15.8 s to +3.0 s.
    ///
    /// The jitter is the point. A *constant* lag is invisible to a constant-velocity
    /// filter — it simply tracks a target that is uniformly a little behind, and the
    /// innovations stay small. It is the variation in age that makes an evenly-moving
    /// aircraft appear to lurch, and that is what a gate is obliged to reject.
    ///
    /// The 0.5 s / 3.5 s alternation is sized to that measured spread rather than
    /// picked for effect. Against this timeline a 200 m/s target produces a worst
    /// innovation of 0.01 sigma when the age is honoured and 16.1 sigma when it is
    /// discarded, so the two tests below sit either side of a wide margin instead of
    /// balancing on the gate.
    ///
    /// Long enough, too, for the velocity covariance to converge: until it does the
    /// gate is several hundred metres wide and waves mis-timed plots through, which is
    /// a fact about a cold filter rather than about timestamps.
    fn jittered_timeline() -> Vec<(u64, f64, f64)> {
        (0..30)
            .map(|step| {
                let arrival = 4 * step + 4;
                let age = if step % 2 == 0 { 0.5 } else { 3.5 };
                (arrival, age, arrival as f64 - age)
            })
            .collect()
    }

    /// 200 m/s northbound from the frame origin at `t` seconds.
    fn northbound_lat(t: f64) -> f64 {
        38.0 + 0.0018 * t
    }

    #[test]
    fn a_jittered_plot_is_placed_at_the_moment_it_describes() {
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let start = Instant::now();
        for (arrival, age, obs_t) in jittered_timeline() {
            store.ingest(
                &[aged("bbb", northbound_lat(obs_t), 23.7, 10_000.0, age)],
                &f,
                start + Duration::from_secs(arrival),
            );
        }
        let t = store.get("bbb").unwrap();
        assert_eq!(
            t.superseded, 0,
            "correctly aged plots must not be superseded"
        );
        assert_eq!(t.rejected, 0, "correctly aged plots must not be gated");
        assert!(
            (t.velocity().n - 200.0).abs() < 25.0,
            "vn was {}",
            t.velocity().n
        );
    }

    #[test]
    fn discarding_the_reported_age_gates_those_same_good_plots() {
        // The regression this change exists for, stated as a test. Same aircraft, same
        // arrivals, same true positions — but every plot claims to be fresh, which is
        // exactly what stamping measurements with the cycle clock amounts to. The
        // filter then watches a steadily-flying aircraft lurch back and forth by 200 m
        // and correctly refuses to believe it.
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let start = Instant::now();
        for (arrival, _age, obs_t) in jittered_timeline() {
            store.ingest(
                &[contact("ccc", northbound_lat(obs_t), 23.7, 10_000.0)],
                &f,
                start + Duration::from_secs(arrival),
            );
        }
        assert!(
            store.get("ccc").unwrap().rejected > 0,
            "mis-timed plots should be gated — if this stops holding, the gate is loose"
        );
    }

    #[test]
    fn a_re_served_snapshot_is_superseded_not_reapplied() {
        // The feed hands back the same underlying observation on the next poll: arrival
        // and reported age have both advanced by 2 s, so it describes the same instant.
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let start = Instant::now();
        store.ingest(
            &[aged("ddd", 38.0, 23.7, 10_000.0, 0.5)],
            &f,
            start + Duration::from_secs(1),
        );
        let report = store.ingest(
            &[aged("ddd", 38.0, 23.7, 10_000.0, 2.5)],
            &f,
            start + Duration::from_secs(3),
        );
        assert_eq!(report.superseded, 1);
        let t = store.get("ddd").unwrap();
        assert_eq!(
            t.updates, 1,
            "a duplicate must not shrink the covariance twice"
        );
        assert_eq!(t.superseded, 1);
        assert_eq!(store.total_superseded, 1);
    }

    #[test]
    fn an_out_of_order_plot_is_superseded_not_retrodicted() {
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let start = Instant::now();
        store.ingest(&[contact("eee", 38.0, 23.7, 10_000.0)], &f, start);
        // Arrives later but describes a moment 10 s before the one already held.
        let report = store.ingest(
            &[aged("eee", 38.02, 23.7, 10_000.0, 12.0)],
            &f,
            start + Duration::from_secs(2),
        );
        assert_eq!(report.superseded, 1);
        assert_eq!(store.get("eee").unwrap().updates, 1);
    }

    #[test]
    fn a_gated_plot_still_coasts_the_track() {
        // Rejecting a measurement is not rejecting the passage of time.
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let start = Instant::now();
        for step in 0..10 {
            store.ingest(
                &[contact("fff", 38.0 + 0.0018 * step as f64, 23.7, 10_000.0)],
                &f,
                start + Duration::from_secs(step),
            );
        }
        let report = store.ingest(
            &[contact("fff", 43.4, 23.7, 10_000.0)],
            &f,
            start + Duration::from_secs(11),
        );
        assert_eq!(report.gated, 1);
        let t = store.get("fff").unwrap();
        assert_eq!(t.rejected, 1);
        assert_eq!(
            t.filter_age_at(start + Duration::from_secs(11)),
            0.0,
            "the gated plot should still have advanced the validity time to its own moment"
        );
    }

    #[test]
    fn staleness_counts_from_observation_not_arrival() {
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let start = Instant::now();
        store.ingest(&[aged("ggg", 38.0, 23.7, 10_000.0, 6.0)], &f, start);
        let t = store.get("ggg").unwrap();
        // Arrived just now, but the position inside it was already six seconds old.
        assert!(
            (t.staleness(start) - 6.0).abs() < 1e-9,
            "{}",
            t.staleness(start)
        );
        // And it therefore ages out of a 45 s window six seconds sooner.
        assert_eq!(
            store.prune(start + Duration::from_secs(40), Duration::from_secs(45)),
            1
        );
    }
}
