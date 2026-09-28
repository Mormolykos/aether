//! Runtime configuration read from a `.env` file next to the binary.
//!
//! No config crate: the parser is fifteen lines, it has no transitive dependencies, and
//! on an edge target every dependency is something you have to justify to a certifier.

use crate::domain::{self, within, ENVELOPE};
use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;
use std::path::Path;

/// Hard ceiling on `MAX_BODY_BYTES`: 64 MiB, eight times the shipped default and some
/// 700 times the measured steady-state response. The body limit exists to bound memory,
/// so the configuration must not be able to raise it far enough to stop bounding it.
pub const MAX_BODY_CEILING: usize = 64 * 1024 * 1024;

#[derive(Debug, Clone)]
pub struct Config {
    pub site_lat: f64,
    pub site_lon: f64,
    pub site_alt_m: f64,
    pub adsb_url: String,
    pub adsb_radius_nm: u32,
    /// Hard ceiling on one response body, in bytes.
    pub max_body_bytes: usize,
    pub poll_ms: u64,
    pub cycle_ms: u64,
    pub track_timeout_s: u64,
    pub horizon_s: f64,
    pub min_horiz_m: f64,
    pub min_vert_m: f64,
    pub process_noise: f64,
    pub meas_var: f64,
    pub gate_sigma: f64,
    pub display_rows: usize,
}

impl Config {
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let raw = std::fs::read_to_string(path)
            .with_context(|| format!("cannot read config at {}", path.display()))?;
        let env = parse(&raw);

        let cfg = Self {
            site_lat: get(&env, "SITE_LAT")?,
            site_lon: get(&env, "SITE_LON")?,
            site_alt_m: get(&env, "SITE_ALT_M")?,
            adsb_url: env
                .get("ADSB_URL")
                .cloned()
                .ok_or_else(|| anyhow!("ADSB_URL missing"))?,
            adsb_radius_nm: get(&env, "ADSB_RADIUS_NM")?,
            max_body_bytes: get(&env, "MAX_BODY_BYTES")?,
            poll_ms: get(&env, "POLL_MS")?,
            cycle_ms: get(&env, "CYCLE_MS")?,
            track_timeout_s: get(&env, "TRACK_TIMEOUT_S")?,
            horizon_s: get(&env, "HORIZON_S")?,
            min_horiz_m: get(&env, "MIN_HORIZ_M")?,
            min_vert_m: get(&env, "MIN_VERT_M")?,
            process_noise: get(&env, "PROCESS_NOISE")?,
            meas_var: get(&env, "MEAS_VAR")?,
            gate_sigma: get(&env, "GATE_SIGMA")?,
            display_rows: get(&env, "DISPLAY_ROWS")?,
        };

        let bad = cfg.violations();
        if !bad.is_empty() {
            return Err(anyhow!(
                "invalid configuration in {}: {}",
                path.display(),
                bad.join("; ")
            ));
        }
        Ok(cfg)
    }

    /// Every value checked against its declared domain, with every violation reported at
    /// once. Parsing only establishes that a value is a number; this establishes that it
    /// is one the process can run on. Empty means the configuration is usable.
    pub fn violations(&self) -> Vec<String> {
        let mut bad = Vec::new();
        let mut check = |ok: bool, key: &str, value: &dyn std::fmt::Display, rule: &str| {
            if !ok {
                bad.push(format!("{key}={value}: {rule}"));
            }
        };

        // The observer is the origin of every transformed position; a bad one poisons
        // all of them.
        let (lat, lon, alt) = (self.site_lat, self.site_lon, self.site_alt_m);
        check(
            within(lat, -90.0, 90.0),
            "SITE_LAT",
            &lat,
            "must be a number from -90 to 90",
        );
        check(
            within(lon, -180.0, 180.0),
            "SITE_LON",
            &lon,
            "must be a number from -180 to 180",
        );
        check(
            within(alt, -ENVELOPE, ENVELOPE),
            "SITE_ALT_M",
            &alt,
            "must be a number of metres within ±1e9",
        );

        // What can be settled here, locally and deterministically, is the form of the
        // request URL: it parses, its scheme is https, and it names a host. Whether that
        // host resolves, answers, or serves the expected JSON is a runtime fact, and it is
        // counted by the sensor-health line, not predicted here. The README's claim that
        // the feed is read over TLS rests on this check plus the client's refusal of any
        // non-https hop.
        let endpoint = self.adsb_endpoint();
        // The parser follows the WHATWG URL rules, which forgive extra slashes: it reads
        // `https:///v2/lat` as host "v2". As written, that URL has an empty authority,
        // which RFC 3986 does not allow for https, so the authority is checked on the
        // text itself as well.
        let authority_written = endpoint
            .get(..8)
            .is_some_and(|s| s.eq_ignore_ascii_case("https://"))
            && endpoint[8..]
                .chars()
                .next()
                .is_some_and(|ch| !matches!(ch, '/' | '\\' | '?' | '#'));
        let usable = authority_written
            && reqwest::Url::parse(&endpoint).is_ok_and(|u| {
                u.scheme() == "https" && u.host_str().is_some_and(|h| !h.is_empty())
            });
        check(
            usable,
            "ADSB_URL",
            &self.adsb_url,
            "must be an absolute https:// URL with a host",
        );
        check(
            self.adsb_radius_nm >= 1,
            "ADSB_RADIUS_NM",
            &self.adsb_radius_nm,
            "must be at least 1",
        );
        // Operating policy, not a derived bound: the ceiling keeps one response's buffer
        // bounded per poll. It says nothing about aggregate memory.
        check(
            (1..=MAX_BODY_CEILING).contains(&self.max_body_bytes),
            "MAX_BODY_BYTES",
            &self.max_body_bytes,
            "must be from 1 to 67108864 (64 MiB, a per-response policy ceiling)",
        );

        // A zero period panics the timer: that part is necessary. Requiring a period
        // shorter than the track timeout is operating policy, not a correctness
        // invariant. A longer period is not wrong in itself; it gives a picture that
        // empties between polls, which this demonstrator treats as misconfiguration.
        let timeout_ms = u128::from(self.track_timeout_s) * 1000;
        check(
            self.track_timeout_s >= 1,
            "TRACK_TIMEOUT_S",
            &self.track_timeout_s,
            "must be at least 1",
        );
        for (key, ms) in [("POLL_MS", self.poll_ms), ("CYCLE_MS", self.cycle_ms)] {
            check(
                ms >= 1 && u128::from(ms) < timeout_ms,
                key,
                &ms,
                "must be at least 1 and, as operating policy, shorter than TRACK_TIMEOUT_S",
            );
        }

        // A NaN minimum makes every comparison false, and a zero one can never be
        // breached: either way the screen goes quiet without saying so.
        for (key, v, why) in [
            ("HORIZON_S", self.horizon_s, domain::horizon(self.horizon_s)),
            (
                "MIN_HORIZ_M",
                self.min_horiz_m,
                domain::separation_minimum(self.min_horiz_m),
            ),
            (
                "MIN_VERT_M",
                self.min_vert_m,
                domain::separation_minimum(self.min_vert_m),
            ),
        ] {
            if let Some(rule) = why {
                check(false, key, &v, rule);
            }
        }

        // The same predicates `TrackStore` enforces, so a configuration cannot pass here
        // and then be refused when the store is built.
        for (key, v, why) in [
            (
                "PROCESS_NOISE",
                self.process_noise,
                domain::process_noise(self.process_noise),
            ),
            ("MEAS_VAR", self.meas_var, domain::meas_var(self.meas_var)),
            (
                "GATE_SIGMA",
                self.gate_sigma,
                domain::gate_sigma(self.gate_sigma),
            ),
        ] {
            if let Some(rule) = why {
                check(false, key, &v, rule);
            }
        }
        bad
    }

    /// Expand the feed template against the configured site.
    pub fn adsb_endpoint(&self) -> String {
        self.adsb_url
            .replace("{lat}", &format!("{:.4}", self.site_lat))
            .replace("{lon}", &format!("{:.4}", self.site_lon))
            .replace("{dist}", &self.adsb_radius_nm.to_string())
    }
}

fn parse(raw: &str) -> HashMap<String, String> {
    let mut out = HashMap::new();
    for line in raw.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let Some((key, value)) = line.split_once('=') else {
            continue;
        };
        let value = value.trim().trim_matches('"').trim_matches('\'');
        out.insert(key.trim().to_string(), value.to_string());
    }
    out
}

fn get<T: std::str::FromStr>(env: &HashMap<String, String>, key: &str) -> Result<T> {
    let raw = env
        .get(key)
        .ok_or_else(|| anyhow!("{key} missing from .env"))?;
    raw.parse::<T>()
        .map_err(|_| anyhow!("{key}={raw} is not a valid {}", std::any::type_name::<T>()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conjunction;
    use crate::geo::{Frame, Geodetic};
    use crate::ingest::Contact;
    use crate::track::TrackStore;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::time::{Duration, Instant};

    const EXAMPLE: &str = include_str!("../.env.example");

    /// The shipped example with some keys replaced, loaded through the real file path.
    fn load_with(overrides: &[(&str, &str)]) -> Result<Config> {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let mut text = String::new();
        for line in EXAMPLE.lines() {
            let key = line.split_once('=').map(|(k, _)| k.trim());
            match overrides.iter().find(|(k, _)| Some(*k) == key) {
                Some((k, v)) => text.push_str(&format!("{k}={v}\n")),
                None => text.push_str(&format!("{line}\n")),
            }
        }
        let path = std::env::temp_dir().join(format!(
            "aether-config-test-{}-{}.env",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, text).expect("write test config");
        let loaded = Config::load(&path);
        let _ = std::fs::remove_file(&path);
        loaded
    }

    #[test]
    fn the_example_configuration_is_accepted() {
        load_with(&[]).expect("the shipped example must load");
    }

    #[test]
    fn configuration_that_would_crash_or_disable_checks_is_refused() {
        // Audit regression: every case the audit named, each alone on an otherwise
        // valid file, and the refusal must name the key at fault.
        for (key, value) in [
            ("POLL_MS", "0"),
            ("CYCLE_MS", "0"),
            ("TRACK_TIMEOUT_S", "0"),
            ("HORIZON_S", "-1"),
            ("HORIZON_S", "NaN"),
            ("HORIZON_S", "inf"),
            ("GATE_SIGMA", "NaN"),
            ("GATE_SIGMA", "0"),
            ("GATE_SIGMA", "-5"),
            ("PROCESS_NOISE", "-1"),
            ("PROCESS_NOISE", "NaN"),
            ("MEAS_VAR", "-1"),
            ("MEAS_VAR", "0"),
            ("MIN_HORIZ_M", "NaN"),
            ("MIN_HORIZ_M", "0"),
            ("MIN_VERT_M", "-305"),
            ("MIN_VERT_M", "inf"),
            ("SITE_LAT", "91"),
            ("SITE_LAT", "NaN"),
            ("SITE_LON", "-180.5"),
            ("SITE_ALT_M", "inf"),
            ("SITE_ALT_M", "1e300"),
            ("MAX_BODY_BYTES", "0"),
            ("MAX_BODY_BYTES", "18446744073709551615"),
            (
                "ADSB_URL",
                "http://api.adsb.lol/v2/lat/{lat}/lon/{lon}/dist/{dist}",
            ),
            // Subnormal values: fewer than 53 significant bits, so the value used is not
            // the value written (3e-324 is read as 5e-324).
            ("MEAS_VAR", "5e-324"),
            ("MIN_HORIZ_M", "3e-324"),
            ("MIN_VERT_M", "1e-310"),
            ("GATE_SIGMA", "1e-310"),
            ("HORIZON_S", "1e-320"),
            ("PROCESS_NOISE", "1e-315"),
            // Audit regression: the scheme alone, with no host, passed a prefix check.
            ("ADSB_URL", "https://"),
            ("ADSB_URL", "https:///v2/lat/{lat}"),
            ("ADSB_URL", "https://exa mple.com/"),
            ("ADSB_URL", "https//api.adsb.lol/"),
            ("ADSB_URL", "ftp://api.adsb.lol/"),
            ("ADSB_RADIUS_NM", "0"),
        ] {
            match load_with(&[(key, value)]) {
                Ok(_) => panic!("{key}={value} was accepted"),
                Err(e) => assert!(
                    e.to_string().contains(key),
                    "{key}={value} refused, but the error does not name it: {e}"
                ),
            }
        }
        // Zero process noise is a legitimate tuning on its own; with zero measurement
        // variance as well, the innovation variance collapses to zero.
        assert!(load_with(&[("PROCESS_NOISE", "0"), ("MEAS_VAR", "0")]).is_err());
    }

    #[test]
    fn legitimate_edge_values_are_accepted() {
        for overrides in [
            &[("HORIZON_S", "0")][..],
            &[("PROCESS_NOISE", "0")][..],
            &[("SITE_LAT", "-90"), ("SITE_LON", "180")][..],
            &[("SITE_ALT_M", "-400")][..],
            &[("MAX_BODY_BYTES", "65536")][..],
            &[("DISPLAY_ROWS", "0")][..],
            // Tiny but normal: the arithmetic is written to be right for these.
            &[("MIN_HORIZ_M", "1e-300"), ("MIN_VERT_M", "1e-300")][..],
            &[("MEAS_VAR", "1e-20"), ("PROCESS_NOISE", "0")][..],
            &[(
                "ADSB_URL",
                "HTTPS://api.adsb.lol/v2/lat/{lat}/lon/{lon}/dist/{dist}",
            )][..],
        ] {
            if let Err(e) = load_with(overrides) {
                panic!("{overrides:?} was refused: {e}");
            }
        }
    }

    #[test]
    fn an_accepted_tiny_minimum_still_catches_colocation() {
        // Audit regression, as the audit ran it: a configuration with a horizontal
        // minimum of 1e-300 m is accepted, and two tracks in exactly the same place must
        // then breach it. Squaring that minimum underflowed to 0, and the screen asked
        // whether 0 < 0.
        let cfg = load_with(&[("MIN_HORIZ_M", "1e-300")]).expect("a normal positive minimum");
        let frame = Frame::new(Geodetic {
            lat_deg: 0.0,
            lon_deg: 0.0,
            alt_m: 0.0,
        });
        let still = |id: &str| Contact {
            id: id.into(),
            label: id.into(),
            kind: String::new(),
            squawk: String::new(),
            emergency: false,
            on_ground: false,
            geo: Geodetic {
                lat_deg: 0.0,
                lon_deg: 0.0,
                alt_m: 0.0,
            },
            gs_kt: None,
            track_deg: None,
            vrate_fpm: None,
            age_s: 0.0,
            source: "TEST",
        };
        let mut store = TrackStore::new(cfg.process_noise, cfg.meas_var, cfg.gate_sigma);
        let now = Instant::now();
        store.ingest(&[still("a"), still("b")], &frame, now);
        let view = |id: &str| store.get(id).unwrap().view_at(now);
        let alert = conjunction::pair_cpa(
            store.get("a").unwrap(),
            view("a"),
            store.get("b").unwrap(),
            view("b"),
            cfg.horizon_s,
            cfg.min_horiz_m,
            cfg.min_vert_m,
        );
        assert!(
            alert.is_some(),
            "colocated tracks inside an accepted minimum"
        );
    }

    #[test]
    fn every_violation_is_reported_at_once() {
        let e = load_with(&[("POLL_MS", "0"), ("GATE_SIGMA", "NaN")])
            .expect_err("two violations")
            .to_string();
        assert!(e.contains("POLL_MS") && e.contains("GATE_SIGMA"), "{e}");
    }

    /// Runs what the process does with a configuration, without the network: builds the
    /// frame and the store, feeds a firm track and a head-on pair, screens them. Returns
    /// what went wrong, if anything. It shares no code with the validator, so it is an
    /// independent witness of whether an accepted configuration is actually usable.
    fn exercise(cfg: &Config) -> Option<String> {
        for (name, ms) in [("poll", cfg.poll_ms), ("cycle", cfg.cycle_ms)] {
            if Duration::from_millis(ms).is_zero() {
                return Some(format!("{name} period is zero"));
            }
        }
        let cfg = cfg.clone();
        std::panic::catch_unwind(move || {
            let frame = Frame::new(Geodetic {
                lat_deg: cfg.site_lat,
                lon_deg: cfg.site_lon,
                alt_m: cfg.site_alt_m,
            });
            let mut store = TrackStore::new(cfg.process_noise, cfg.meas_var, cfg.gate_sigma);
            let start = Instant::now();
            let plot = |id: &str, lat: f64| Contact {
                id: id.into(),
                label: id.into(),
                kind: String::new(),
                squawk: String::new(),
                emergency: false,
                on_ground: false,
                geo: Geodetic {
                    lat_deg: lat,
                    lon_deg: cfg.site_lon,
                    alt_m: 10_000.0,
                },
                gs_kt: None,
                track_deg: None,
                vrate_fpm: None,
                age_s: 0.0,
                source: "TEST",
            };
            for step in 0..12u64 {
                let s = step as f64;
                store.ingest(
                    &[
                        plot("north", cfg.site_lat + 0.2 - 0.0018 * s),
                        plot("south", cfg.site_lat - 0.2 + 0.0018 * s),
                    ],
                    &frame,
                    start + Duration::from_secs(step),
                );
            }
            let at = start + Duration::from_secs(11);
            for t in store.iter() {
                let v = t.view_at(at);
                let all = [v.pos.e, v.pos.n, v.pos.u, v.vel.e, v.vel.n, v.vel.u];
                if all.iter().any(|x| !x.is_finite()) {
                    return Some(format!("track {} went non-finite", t.id));
                }
            }
            let alerts =
                conjunction::screen(&store, at, cfg.horizon_s, cfg.min_horiz_m, cfg.min_vert_m);
            for c in &alerts {
                let all = [c.t_los, c.t_cpa, c.horiz_m, c.vert_m, c.closing_ms];
                if all.iter().any(|x| !x.is_finite()) {
                    return Some("an alert carried a non-finite number".into());
                }
            }
            None
        })
        .unwrap_or_else(|_| Some("panicked".into()))
    }

    #[test]
    fn no_single_bad_value_is_both_accepted_and_unusable() {
        // The property rather than the list: every numeric key, set to every value a
        // hand-edited file can plausibly contain, is either refused by name or yields a
        // configuration the process can actually run.
        let keys = [
            "SITE_LAT",
            "SITE_LON",
            "SITE_ALT_M",
            "ADSB_RADIUS_NM",
            "MAX_BODY_BYTES",
            "POLL_MS",
            "CYCLE_MS",
            "TRACK_TIMEOUT_S",
            "HORIZON_S",
            "MIN_HORIZ_M",
            "MIN_VERT_M",
            "PROCESS_NOISE",
            "MEAS_VAR",
            "GATE_SIGMA",
            "DISPLAY_ROWS",
        ];
        let values = [
            "NaN",
            "inf",
            "-inf",
            "0",
            "-0",
            "-1",
            "1",
            "0.5",
            "1e-300",
            "1e300",
            "-1e300",
            "1e9",
            "180",
            "-180",
            "90.0001",
            "18446744073709551615",
        ];
        let mut failures = Vec::new();
        for key in keys {
            for value in values {
                match load_with(&[(key, value)]) {
                    Err(e) if !e.to_string().contains(key) => {
                        failures.push(format!("{key}={value}: refused without naming it: {e}"))
                    }
                    Err(_) => {}
                    Ok(cfg) => {
                        if let Some(why) = exercise(&cfg) {
                            failures.push(format!("{key}={value}: accepted, then {why}"));
                        }
                    }
                }
            }
        }
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    }
}
