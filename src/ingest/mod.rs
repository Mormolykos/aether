//! Sensor adapters.
//!
//! Every adapter converts whatever wire format it speaks into `Contact` and pushes
//! batches down an mpsc channel. Nothing downstream knows whether a contact arrived
//! from an ADS-B feed, a radar extractor or a TLE catalogue, so adding a sensor is
//! adding one file here and one spawn in `main`.

pub mod adsb;

use crate::geo::Geodetic;
use std::time::Instant;

/// One poll's worth of contacts, carrying the instant the data actually landed.
///
/// The arrival instant belongs to the batch rather than to each contact because it is a
/// property of the transfer, not of the observation. Each contact then says how old it
/// already was at that moment (`age_s`), and the tracker works entirely in differences
/// of the two — see `Track::filter_age_at`.
///
/// A failed poll still produces a batch, with no contacts and an updated `health`. That
/// is the whole reason health rides here rather than in a shared counter: it keeps the
/// channel the single path out of the sensor task, so there is no `Arc`, no lock, and
/// no second thing to keep in sync. An empty batch is not an observation of an empty
/// sky — `health` is what distinguishes the two.
#[derive(Debug)]
pub struct Batch {
    pub received: Instant,
    pub contacts: Vec<Contact>,
    pub health: SensorHealth,
}

/// Cumulative outcome counters for one sensor adapter.
///
/// Owned by the sensor task and copied onto every batch it emits. Small and `Copy` on
/// purpose: sending a snapshot costs nothing and cannot go stale in a way that matters.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SensorHealth {
    /// Polls attempted.
    pub polls: u64,
    /// Transport failures and non-success statuses.
    pub http_errors: u64,
    /// Responses refused for exceeding the configured body limit.
    pub oversized_responses: u64,
    /// Bodies that arrived within the limit but were not the JSON we expect.
    pub decode_errors: u64,
}

impl SensorHealth {
    /// Total failed polls.
    pub fn faults(&self) -> u64 {
        self.http_errors + self.oversized_responses + self.decode_errors
    }

    fn record(&mut self, err: &IngestError) {
        match err {
            IngestError::Http(_) => self.http_errors += 1,
            IngestError::Oversized { .. } => self.oversized_responses += 1,
            IngestError::Decode(_) => self.decode_errors += 1,
        }
    }
}

/// Why one poll produced nothing usable.
///
/// Typed rather than a string so the sensor can count the kinds separately: a feed that
/// is unreachable, a feed that is flooding us, and a feed that has changed its schema
/// are three different problems with three different responses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IngestError {
    /// Transport failure, or a non-success status from the feed.
    Http(String),
    /// The response exceeded the configured body limit and was abandoned.
    Oversized { limit: usize, seen: usize },
    /// The body arrived intact but did not parse.
    Decode(String),
}

impl std::fmt::Display for IngestError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Http(e) => write!(f, "http: {e}"),
            Self::Oversized { limit, seen } => {
                write!(f, "response body exceeded {limit} bytes (reached {seen})")
            }
            Self::Decode(e) => write!(f, "decode: {e}"),
        }
    }
}

impl std::error::Error for IngestError {}

/// One plot from one sensor at one instant, before association.
#[derive(Clone, Debug)]
pub struct Contact {
    /// Sensor-scoped identity used for association. ICAO 24-bit address for ADS-B.
    pub id: String,
    pub label: String,
    pub kind: String,
    pub squawk: String,
    pub emergency: bool,
    pub on_ground: bool,
    pub geo: Geodetic,
    /// Reported ground speed, knots. Reported, not estimated: the filter owns velocity.
    pub gs_kt: Option<f64>,
    pub track_deg: Option<f64>,
    pub vrate_fpm: Option<f64>,
    /// How many seconds old this position already was when the batch arrived, as
    /// reported by the sensor itself. Not a derived quantity and not an arrival time:
    /// keeping it relative is what lets the tracker avoid reconstructing an absolute
    /// observation instant.
    pub age_s: f64,
    pub source: &'static str,
}
