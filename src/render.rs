//! Operator picture and cycle telemetry.
//!
//! The whole frame is composed into one buffer and written with a single call. Printing
//! field by field is how a display task ends up owning the frame time of the process
//! it is supposed to be observing.

use crate::conjunction::Conjunction;
use crate::geo::Frame;
use crate::ingest::SensorHealth;
use crate::track::TrackStore;
use std::fmt::Write as _;
use std::io::Write as _;
use std::time::{Duration, Instant};

/// Rolling cycle telemetry. A system that cannot report its own timing cannot be
/// trusted to hold a deadline.
pub struct Cycles {
    pub started: Instant,
    pub count: u64,
    pub last: Duration,
    pub worst: Duration,
    total: Duration,
    pub contacts_in: u64,
    pub polls_seen: u64,
}

impl Default for Cycles {
    fn default() -> Self {
        Self {
            started: Instant::now(),
            count: 0,
            last: Duration::ZERO,
            worst: Duration::ZERO,
            total: Duration::ZERO,
            contacts_in: 0,
            polls_seen: 0,
        }
    }
}

impl Cycles {
    pub fn record(&mut self, elapsed: Duration, contacts: usize, polls: usize) {
        self.count += 1;
        self.last = elapsed;
        self.worst = self.worst.max(elapsed);
        self.total += elapsed;
        self.contacts_in += contacts as u64;
        self.polls_seen += polls as u64;
    }

    pub fn mean(&self) -> Duration {
        if self.count == 0 {
            Duration::ZERO
        } else {
            self.total / self.count as u32
        }
    }
}

const CLEAR: &str = "\x1b[2J\x1b[H";
const DIM: &str = "\x1b[2m";
const BOLD: &str = "\x1b[1m";
const RED: &str = "\x1b[31m";
const YELLOW: &str = "\x1b[33m";
const CYAN: &str = "\x1b[36m";
const RESET: &str = "\x1b[0m";

pub fn draw(
    store: &TrackStore,
    alerts: &[Conjunction],
    cyc: &Cycles,
    health: &SensorHealth,
    frame: &Frame,
    rows: usize,
    now: Instant,
) {
    let mut s = String::with_capacity(8192);
    s.push_str(CLEAR);

    let up = now.saturating_duration_since(cyc.started).as_secs();
    let _ = writeln!(
        s,
        "{BOLD}AETHER{RESET}  air picture   up {:02}:{:02}:{:02}   cycle {}",
        up / 3600,
        (up / 60) % 60,
        up % 60,
        cyc.count
    );
    let _ = writeln!(
        s,
        "{DIM}tracks {:<4} initiated {:<5} dropped {:<5} gated {:<4} superseded {:<5} \
         plots {:<7} cycle last {:.2}ms mean {:.2}ms worst {:.2}ms{RESET}",
        store.len(),
        store.total_initiated,
        store.total_dropped,
        store.total_gated,
        store.total_superseded,
        cyc.contacts_in,
        cyc.last.as_secs_f64() * 1000.0,
        cyc.mean().as_secs_f64() * 1000.0,
        cyc.worst.as_secs_f64() * 1000.0,
    );

    // Sensor health is stated separately from the tracker's counters, because an empty
    // picture caused by an empty sky and an empty picture caused by a sensor that has
    // stopped answering look identical everywhere else on this screen.
    let sensor_colour = if health.faults() > 0 { YELLOW } else { DIM };
    let _ = writeln!(
        s,
        "{sensor_colour}sensor  polls {:<6} http-err {:<4} oversized {:<4} decode-err {:<4}{RESET}",
        health.polls, health.http_errors, health.oversized_responses, health.decode_errors,
    );
    s.push('\n');

    let _ = writeln!(
        s,
        "{BOLD}{:<9} {:<9} {:<6} {:>7} {:>6} {:>7} {:>5} {:>7} {:>6} {:>6} {:>5}{RESET}",
        "ID", "CALLSIGN", "TYPE", "RNG km", "BRG", "ALT ft", "GS kt", "HDG", "V fpm", "±m", "AGE"
    );

    // Every track extrapolated to the same instant. These are views: drawing the
    // picture must not advance anybody's filter.
    let mut sorted: Vec<_> = store.iter().map(|t| (t, t.view_at(now))).collect();
    sorted.sort_by(|(_, a), (_, b)| a.pos.horiz().total_cmp(&b.pos.horiz()));

    for (t, v) in sorted.iter().take(rows) {
        let p = v.pos;
        let flag = if t.emergency {
            RED
        } else if !t.is_firm() {
            DIM
        } else {
            ""
        };
        let _ = writeln!(
            s,
            "{flag}{:<9} {:<9} {:<6} {:>7.1} {:>5.0}° {:>7.0} {:>5.0} {:>6.0}° {:>6.0} {:>6.0} {:>4.0}s{RESET}",
            t.id,
            truncate(&t.label, 9),
            truncate(&t.kind, 6),
            p.horiz() / 1000.0,
            p.bearing_deg(),
            // True altitude above the ellipsoid, not the ENU "Up" axis. The two agree
            // near the origin and diverge as the square of the range.
            frame.to_geodetic(p).alt_m / crate::geo::FT_TO_M,
            t.speed_kt(),
            t.heading_deg(),
            t.vrate_fpm(),
            v.pos_sigma,
            t.staleness(now),
        );
    }
    if store.len() > rows {
        let _ = writeln!(s, "{DIM}   ... {} more{RESET}", store.len() - rows);
    }

    s.push('\n');
    if alerts.is_empty() {
        let _ = writeln!(s, "{DIM}conjunction screen clear{RESET}");
    } else {
        let _ = writeln!(
            s,
            "{BOLD}{YELLOW}CONJUNCTION SCREEN — {} pair(s){RESET}",
            alerts.len()
        );
        for c in alerts.iter().take(6) {
            let colour = if c.t_cpa < 60.0 { RED } else { YELLOW };
            let _ = writeln!(
                s,
                "{colour}  T-{:>3.0}s  {:<9} <-> {:<9}  horiz {:>5.2} km  vert {:>4.0} m  \
                 closing {:>3.0} kt  ±{:.0} m{RESET}",
                c.t_cpa,
                truncate(&c.label_a, 9),
                truncate(&c.label_b, 9),
                c.horiz_m / 1000.0,
                c.vert_m,
                c.closing_ms * crate::geo::MS_TO_KT,
                c.sigma_m,
            );
        }
    }

    let _ = writeln!(s, "\n{CYAN}{DIM}ctrl-c to stop{RESET}");

    let stdout = std::io::stdout();
    let mut lock = stdout.lock();
    let _ = lock.write_all(s.as_bytes());
    let _ = lock.flush();
}

fn truncate(s: &str, n: usize) -> &str {
    match s.char_indices().nth(n) {
        Some((idx, _)) => &s[..idx],
        None => s,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn truncate_is_char_safe() {
        assert_eq!(truncate("ABCDEFGHIJ", 4), "ABCD");
        assert_eq!(truncate("AB", 4), "AB");
        // Must not split a multi-byte character and panic on a slice boundary.
        assert_eq!(truncate("ΑΘΗΝΑ", 3), "ΑΘΗ");
    }

    #[test]
    fn mean_of_no_cycles_is_zero() {
        assert_eq!(Cycles::default().mean(), Duration::ZERO);
    }

    #[test]
    fn worst_case_is_retained() {
        let mut c = Cycles::default();
        c.record(Duration::from_millis(1), 10, 1);
        c.record(Duration::from_millis(9), 10, 1);
        c.record(Duration::from_millis(2), 10, 1);
        assert_eq!(c.worst, Duration::from_millis(9));
        assert_eq!(c.contacts_in, 30);
    }
}
