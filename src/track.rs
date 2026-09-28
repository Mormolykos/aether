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

use crate::domain::{self, ENVELOPE};
use crate::geo::{Enu, Frame, Geodetic};
use crate::ingest::Contact;
use crate::kalman::{Innovation, Kf1D};
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
    /// Outside the declared domain, or would have left the estimate non-finite. The
    /// track, if there is one, is left exactly as it was.
    Invalid,
}

/// Whether a filter state may be kept: every term finite, and no velocity faster than
/// the envelope. A state that fails this came from arithmetic, not from an aircraft, and
/// keeping it would let the next extrapolation overflow.
fn usable(filt: &[Kf1D; 3]) -> bool {
    filt.iter().all(|f| f.is_finite() && f.v.abs() <= ENVELOPE)
}

impl Track {
    /// A new track on its first plot, or `None` if the state it would start from is not
    /// usable.
    fn new(c: &Contact, pos: Enu, received: Instant, q: f64, r: f64) -> Option<Self> {
        // Seed velocity from the reported heading and ground speed when present. It is
        // only a prior: the covariance is loose enough that two updates overrule it. A
        // reported value outside its domain is a value that was not reported.
        let (ve, vn) = match (
            domain::ground_speed_ms(c.gs_kt),
            domain::bearing_rad(c.track_deg),
        ) {
            (Some(speed), Some(rad)) => (speed * rad.sin(), speed * rad.cos()),
            _ => (0.0, 0.0),
        };
        let vu = domain::climb_ms(c.vrate_fpm).unwrap_or(0.0);

        let track = Self {
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
        };
        usable(&track.filt).then_some(track)
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

        // Everything below works on copies, and the track takes a copy only once it is
        // known to be usable. A plot that would leave the estimate non-finite therefore
        // leaves it exactly as it was, as if the plot had never arrived.
        let mut coasted = self.filt;
        for f in &mut coasted {
            f.predict(dt);
        }
        let innovs = [
            coasted[0].innovation(pos.e),
            coasted[1].innovation(pos.n),
            coasted[2].innovation(pos.u),
        ];
        // Checked before the gate and whatever the track's maturity. The gate asks
        // whether a residual is too large, and a NaN is never too large: it has to be
        // refused for being unusable, not left to fail a comparison it cannot fail.
        if !usable(&coasted) || !innovs.iter().all(Innovation::is_usable) {
            return Outcome::Invalid;
        }

        // Coasting to the measurement's own moment is correct whether or not the
        // measurement survives the gate, so the validity time moves first and stays
        // moved. A gated plot leaves a coasted track, which is what it should leave.
        //
        // Gate on all three axes before touching any of them: a partially applied
        // update on a bad plot is worse than a rejected one. Written so that only a
        // sigma inside the gate passes, rather than so that one outside it fails.
        if self.is_firm() && !innovs.iter().all(|i| i.sigma() <= gate_sigma) {
            self.filt = coasted;
            self.valid_rx = received;
            self.valid_age = c.age_s;
            self.rejected += 1;
            return Outcome::Gated;
        }

        let mut updated = coasted;
        for (f, z) in updated.iter_mut().zip([pos.e, pos.n, pos.u]) {
            f.update(z);
        }
        let nis = innovs.iter().map(|i| i.nis()).sum::<f64>() / 3.0;
        if !usable(&updated) || !nis.is_finite() {
            return Outcome::Invalid;
        }

        self.filt = updated;
        self.valid_rx = received;
        self.valid_age = c.age_s;
        self.last_nis = nis;
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
    /// Refused as outside the declared domain or numerically unusable.
    pub invalid: usize,
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
    pub total_invalid: u64,
}

impl TrackStore {
    /// A store with the given filter tuning, or the reasons that tuning is unusable.
    pub fn try_new(process_noise: f64, meas_var: f64, gate_sigma: f64) -> Result<Self, String> {
        let bad: Vec<String> = [
            (
                "process noise",
                process_noise,
                domain::process_noise(process_noise),
            ),
            ("measurement variance", meas_var, domain::meas_var(meas_var)),
            ("gate sigma", gate_sigma, domain::gate_sigma(gate_sigma)),
        ]
        .into_iter()
        .filter_map(|(name, v, why)| why.map(|w| format!("{name} {v}: {w}")))
        .collect();
        if !bad.is_empty() {
            return Err(bad.join("; "));
        }
        Ok(Self {
            tracks: HashMap::new(),
            q: process_noise,
            r: meas_var,
            gate_sigma,
            total_initiated: 0,
            total_dropped: 0,
            total_gated: 0,
            total_superseded: 0,
            total_invalid: 0,
        })
    }

    /// As [`TrackStore::try_new`].
    ///
    /// # Panics
    /// If the tuning is unusable. The configuration is validated before a store is ever
    /// built, so in the application this is unreachable; for a library caller it refuses
    /// at construction rather than producing a store whose gate is silently disabled.
    pub fn new(process_noise: f64, meas_var: f64, gate_sigma: f64) -> Self {
        Self::try_new(process_noise, meas_var, gate_sigma)
            .unwrap_or_else(|e| panic!("unusable filter tuning: {e}"))
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
        // A position is the product of a contact and a frame, so both are checked here,
        // at the only door into the estimator, rather than trusted from whoever built
        // them: `Contact`, `Frame` and this method are all public.
        let frame_ok = frame.is_valid();
        for c in contacts {
            let outcome = if !frame_ok || domain::contact(c).is_some() {
                Outcome::Invalid
            } else {
                let pos = frame.to_enu(c.geo);
                match self.tracks.get_mut(&c.id) {
                    Some(track) => track.update(c, pos, received, self.gate_sigma),
                    None => match Track::new(c, pos, received, self.q, self.r) {
                        Some(track) => {
                            self.tracks.insert(c.id.clone(), track);
                            Outcome::Initiated
                        }
                        None => Outcome::Invalid,
                    },
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
                Outcome::Invalid => {
                    report.invalid += 1;
                    self.total_invalid += 1;
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

    // --- numerical domain (audit finding 7) ------------------------------------------

    /// Every number a track holds, and every number its view hands the screen, is finite.
    fn finite_at(t: &Track, at: Instant) -> bool {
        let v = t.view_at(at);
        let fin = |p: Enu| p.e.is_finite() && p.n.is_finite() && p.u.is_finite();
        fin(t.position())
            && fin(t.velocity())
            && fin(v.pos)
            && fin(v.vel)
            && v.pos_sigma.is_finite()
            && t.pos_sigma().is_finite()
            && t.last_nis.is_finite()
    }

    #[test]
    fn a_huge_finite_ground_speed_cannot_make_the_estimator_non_finite() {
        // Audit regression. A finite ground speed near 1e308 seeds a finite velocity,
        // and multiplying that by a few seconds overflows.
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let start = Instant::now();
        let fast = |lat: f64| Contact {
            gs_kt: Some(1e308),
            ..contact("hhh", lat, 23.7, 10_000.0)
        };
        store.ingest(&[fast(38.0)], &f, start);
        store.ingest(&[fast(38.001)], &f, start + Duration::from_secs(4));
        for t in store.iter() {
            assert!(finite_at(t, start + Duration::from_secs(10)), "{t:?}");
        }
    }

    #[test]
    fn a_non_finite_plot_is_refused_not_absorbed() {
        // Audit regression. The gate asked `sigma > gate`, and NaN compares false, so a
        // NaN innovation passed it. The gate also only runs on firm tracks, so a young
        // track took NaN without being asked at all.
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let start = Instant::now();
        for step in 0..10 {
            store.ingest(
                &[contact("firm", 38.0, 23.7, 10_000.0)],
                &f,
                start + Duration::from_secs(step),
            );
        }
        // One plot: a track too young to be gated at all.
        store.ingest(
            &[contact("new", 38.2, 23.7, 10_000.0)],
            &f,
            start + Duration::from_secs(9),
        );
        let held: Vec<(Enu, u64)> = ["firm", "new"]
            .iter()
            .map(|id| {
                (
                    store.get(id).unwrap().position(),
                    store.get(id).unwrap().updates,
                )
            })
            .collect();

        let at = start + Duration::from_secs(11);
        store.ingest(
            &[
                contact("firm", 38.0, 23.7, f64::NAN),
                contact("new", 38.2, 23.7, f64::NAN),
            ],
            &f,
            at,
        );
        for (id, (pos, updates)) in ["firm", "new"].iter().zip(held) {
            let t = store.get(id).unwrap();
            assert!(finite_at(t, at), "{id}: {t:?}");
            assert_eq!(t.position(), pos, "{id}: the plot moved the estimate");
            assert_eq!(t.updates, updates, "{id}: the plot counted as evidence");
        }
    }

    /// A deterministic generator, so a failing case can be named and replayed.
    struct Lcg(u64);

    impl Lcg {
        fn unit(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
        fn range(&mut self, lo: f64, hi: f64) -> f64 {
            lo + (hi - lo) * self.unit()
        }
        fn pick<T: Copy>(&mut self, xs: &[T]) -> T {
            xs[((self.unit() * xs.len() as f64) as usize).min(xs.len() - 1)]
        }
    }

    /// Values a hostile or broken producer can put in any numeric field.
    const HOSTILE: [f64; 16] = [
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::MAX,
        -f64::MAX,
        1e308,
        -1e308,
        1e200,
        -1e200,
        1e20,
        -1e20,
        1e-300,
        f64::MIN_POSITIVE,
        5e-324,
        -0.0,
        -1.0,
    ];

    #[test]
    fn no_contact_sequence_can_make_a_track_non_finite() {
        // The property rather than the cases: any mix of ordinary and hostile fields,
        // arriving in any order, through the public ingest path that bypasses the
        // decoder. After every batch, every track and every view the screen could take
        // of it, now or long after, is finite.
        let f = frame();
        let start = Instant::now();
        let mut rng = Lcg(0xae7e);
        let (mut refused, mut initiated) = (0, 0);
        for run in 0..40 {
            let mut store = TrackStore::new(4.0, 900.0, 5.0);
            let mut at = start;
            for step in 0..60 {
                at += Duration::from_millis((rng.range(0.0, 3000.0)) as u64);
                let batch: Vec<Contact> = (0..4)
                    .map(|k| {
                        let mut c = Contact {
                            gs_kt: Some(rng.range(0.0, 600.0)),
                            track_deg: Some(rng.range(0.0, 360.0)),
                            vrate_fpm: Some(rng.range(-3000.0, 3000.0)),
                            age_s: rng.range(0.0, 10.0),
                            ..contact(
                                ["p", "q", "r", "s"][k],
                                rng.range(37.5, 38.5),
                                rng.range(23.0, 24.5),
                                rng.range(0.0, 12_000.0),
                            )
                        };
                        // Roughly one contact in four carries at least one bad field.
                        if rng.unit() < 0.25 {
                            let bad = rng.pick(&HOSTILE);
                            match (rng.unit() * 7.0) as u32 {
                                0 => c.geo.lat_deg = bad,
                                1 => c.geo.lon_deg = bad,
                                2 => c.geo.alt_m = bad,
                                3 => c.gs_kt = Some(bad),
                                4 => c.track_deg = Some(bad),
                                5 => c.vrate_fpm = Some(bad),
                                _ => c.age_s = bad,
                            }
                        }
                        c
                    })
                    .collect();
                store.ingest(&batch, &f, at);

                for t in store.iter() {
                    for later in [0, 60, 3_600, 1_000_000] {
                        let when = at + Duration::from_secs(later);
                        assert!(
                            finite_at(t, when),
                            "run {run} step {step}, viewed {later} s on: {t:?}"
                        );
                    }
                    let v = t.velocity();
                    assert!(
                        [v.e, v.n, v.u].iter().all(|x| x.abs() <= ENVELOPE),
                        "run {run} step {step}: velocity beyond the envelope: {t:?}"
                    );
                }
            }
            refused += store.total_invalid;
            initiated += store.total_initiated;
        }
        // Neither half may be vacuous: hostile plots were actually refused, and
        // ordinary ones actually made tracks.
        assert!(
            refused > 100 && initiated > 100,
            "refused {refused}, initiated {initiated}"
        );
    }

    #[test]
    fn stationary_reports_stay_usable_after_a_very_precise_one() {
        // Audit regression, the exact case, through the store. The subtractive covariance
        // update committed a negative velocity variance, after which every identical
        // stationary report was refused as Invalid and the view showed zero uncertainty.
        let origin = Geodetic {
            lat_deg: 0.0,
            lon_deg: 0.0,
            alt_m: 0.0,
        };
        let f = Frame::new(origin);
        let still = Contact {
            gs_kt: None,
            track_deg: None,
            vrate_fpm: None,
            ..contact("still", 0.0, 0.0, 0.0)
        };
        let mut store = TrackStore::new(0.0, 1e-20, 5.0);
        let at = Instant::now();
        store.ingest(std::slice::from_ref(&still), &f, at);
        assert_eq!(
            store
                .ingest(
                    std::slice::from_ref(&still),
                    &f,
                    at + Duration::from_millis(9)
                )
                .updated,
            1
        );
        for s in 1..=5 {
            let when = at + Duration::from_secs(s);
            let report = store.ingest(std::slice::from_ref(&still), &f, when);
            assert_eq!(report.updated, 1, "second {s}: {report:?}");
            let v = store
                .get("still")
                .unwrap()
                .view_at(when + Duration::from_secs(1));
            assert!(
                v.pos_sigma > 0.0,
                "second {s}: uncertainty shown as zero: {v:?}"
            );
        }
    }

    #[test]
    fn a_finite_frame_outside_the_domain_is_refused() {
        // Audit regression. Every field of this frame is finite, and its origin is 1e200 m
        // above the ellipsoid; a valid contact at the ellipsoid then stored an Up
        // coordinate of -1e200.
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let at = Instant::now();
        for site in [
            Geodetic {
                lat_deg: 0.0,
                lon_deg: 0.0,
                alt_m: 1e200,
            },
            Geodetic {
                lat_deg: 1e10,
                lon_deg: 0.0,
                alt_m: 0.0,
            },
            Geodetic {
                lat_deg: 0.0,
                lon_deg: -540.0,
                alt_m: 0.0,
            },
        ] {
            let report = store.ingest(&[contact("a", 0.0, 0.0, 0.0)], &Frame::new(site), at);
            assert_eq!(report.invalid, 1, "{site:?}: {report:?}");
        }
        assert!(store.is_empty());
    }

    #[test]
    fn a_contact_outside_the_domain_is_refused_and_counted() {
        // Direct callers of `ingest`, which never pass through a decoder.
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let f = frame();
        let now = Instant::now();
        let report = store.ingest(
            &[
                contact("lat", 1e300, 23.7, 10_000.0),
                contact("alt", 38.0, 23.7, 1e12),
                aged("age", 38.0, 23.7, 10_000.0, -5.0),
                aged("nan", 38.0, 23.7, 10_000.0, f64::NAN),
            ],
            &f,
            now,
        );
        assert_eq!(report.invalid, 4);
        assert_eq!(report.initiated, 0);
        assert!(store.is_empty());
        assert_eq!(store.total_invalid, 4);
    }

    #[test]
    fn a_broken_frame_cannot_reach_the_filter() {
        // `Frame` is public, and a frame built on a NaN site turns every valid contact
        // into a NaN position after the contact check has passed. The frame check refuses
        // it at the door, on young and firm tracks alike, and the refused plot must leave
        // the track exactly as it was.
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let good = frame();
        let broken = Frame::new(Geodetic {
            lat_deg: f64::NAN,
            lon_deg: 23.7275,
            alt_m: 0.0,
        });
        let start = Instant::now();
        for step in 0..10 {
            store.ingest(
                &[contact("firm", 38.0, 23.7, 10_000.0)],
                &good,
                start + Duration::from_secs(step),
            );
        }
        store.ingest(&[contact("young", 38.1, 23.7, 10_000.0)], &good, start);

        let before: Vec<Track> = ["firm", "young"]
            .iter()
            .map(|id| store.get(id).unwrap().clone())
            .collect();
        let at = start + Duration::from_secs(11);
        let report = store.ingest(
            &[
                contact("firm", 38.0, 23.7, 10_000.0),
                contact("young", 38.1, 23.7, 10_000.0),
                contact("fresh", 38.2, 23.7, 10_000.0),
            ],
            &broken,
            at,
        );
        assert_eq!(report.invalid, 3, "{report:?}");
        assert!(store.get("fresh").is_none());
        for old in before {
            let now = store.get(&old.id).unwrap();
            assert_eq!(now.position(), old.position(), "{}", old.id);
            assert_eq!(now.rejected, old.rejected, "{}: counted as gated", old.id);
            assert_eq!(
                now.filter_age_at(at),
                old.filter_age_at(at),
                "{}: coasted",
                old.id
            );
        }
    }

    #[test]
    fn a_non_finite_position_that_reaches_the_filter_is_refused() {
        // The innovation check, tested directly. With contacts and frames both checked at
        // the door, no public path hands `update` a non-finite position today; this guard
        // is what holds if one ever does, and it must hold on young and firm tracks alike
        // rather than relying on the gate, which asks a question NaN always answers "no" to.
        let f = frame();
        let at = Instant::now();
        let c = contact("t", 38.0, 23.7, 10_000.0);
        let mut young = Track::new(&c, f.to_enu(c.geo), at, 4.0, 900.0).unwrap();
        let mut firm = young.clone();
        firm.updates = 10;
        for track in [&mut young, &mut firm] {
            let before = format!("{track:?}");
            for bad in [f64::NAN, f64::INFINITY] {
                let pos = Enu {
                    e: bad,
                    n: 0.0,
                    u: 0.0,
                };
                let outcome = track.update(&c, pos, at + Duration::from_secs(1), 5.0);
                assert_eq!(outcome, Outcome::Invalid, "{bad}");
                assert_eq!(format!("{track:?}"), before, "{bad}: the track changed");
            }
        }
    }

    #[test]
    fn a_velocity_no_object_can_have_is_never_kept() {
        // A measurement variance of 1e-300 m² is inside the tuning domain, and a time
        // step of 1e-150 s is what two ages 1e-150 s apart produce. Together they give
        // a velocity gain near 1e150, and a 10 km residual becomes 1e154 m/s: finite, so
        // no finiteness check refuses it, and enough to overflow the first long view.
        let mut store = TrackStore::new(4.0, 1e-300, 5.0);
        let f = frame();
        let at = Instant::now();
        store.ingest(&[aged("vvv", 38.0, 23.7, 10_000.0, 1e-150)], &f, at);
        let report = store.ingest(&[aged("vvv", 38.09, 23.7, 10_000.0, 0.0)], &f, at);
        assert_eq!(report.invalid, 1, "{report:?}");
        let v = store.get("vvv").unwrap().velocity();
        assert!([v.e, v.n, v.u].iter().all(|x| x.abs() <= ENVELOPE), "{v:?}");
    }

    #[test]
    fn an_unusable_optional_field_is_absent_not_fatal() {
        // A report with a nonsense ground speed still has a good position: the speed is
        // dropped, the track starts from a zero velocity prior, the position is kept.
        let mut store = TrackStore::new(4.0, 900.0, 5.0);
        let now = Instant::now();
        let report = store.ingest(
            &[Contact {
                gs_kt: Some(1e308),
                vrate_fpm: Some(f64::INFINITY),
                ..contact("opt", 38.0, 23.7, 10_000.0)
            }],
            &frame(),
            now,
        );
        assert_eq!(report.initiated, 1);
        assert_eq!(store.get("opt").unwrap().velocity(), Enu::default());
    }

    #[test]
    fn unusable_tuning_is_refused_at_construction() {
        for (q, r, gate) in [
            (f64::NAN, 900.0, 5.0),
            (-1.0, 900.0, 5.0),
            (4.0, 0.0, 5.0),
            (0.0, 0.0, 5.0),
            (4.0, -1.0, 5.0),
            (4.0, 900.0, f64::NAN),
            (4.0, 900.0, 0.0),
            (4.0, 900.0, f64::INFINITY),
            // Subnormal: not the value written, and a variance that small has no
            // representable covariance.
            (4.0, 5e-324, 5.0),
            (1e-310, 900.0, 5.0),
            (4.0, 900.0, 1e-320),
        ] {
            assert!(TrackStore::try_new(q, r, gate).is_err(), "{q} {r} {gate}");
        }
        assert!(
            TrackStore::try_new(0.0, 900.0, 5.0).is_ok(),
            "zero process noise is legitimate"
        );
        assert!(
            TrackStore::try_new(f64::MIN_POSITIVE, f64::MIN_POSITIVE, f64::MIN_POSITIVE).is_ok(),
            "the smallest normal double is inside the domain"
        );
    }

    #[test]
    #[should_panic(expected = "unusable filter tuning")]
    fn new_panics_rather_than_building_a_store_with_a_disabled_gate() {
        let _ = TrackStore::new(4.0, 900.0, f64::NAN);
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
