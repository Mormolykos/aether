//! Runtime configuration read from a `.env` file next to the binary.
//!
//! No config crate: the parser is fifteen lines, it has no transitive dependencies, and
//! on an edge target every dependency is something you have to justify to a certifier.

use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;
use std::path::Path;

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

        Ok(Self {
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
        })
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
