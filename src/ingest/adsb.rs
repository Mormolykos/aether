//! ADS-B sensor adapter.
//!
//! ADS-B is a cooperative surveillance broadcast: every transponder-equipped aircraft
//! reports its own identity, position, altitude and velocity roughly once a second.
//! It is the cleanest public stand-in for a real air-picture feed, and it has all the
//! properties that make live telemetry awkward: fields go missing, altitude is
//! sometimes the string "ground", positions are stale by an unknown amount, and the
//! producer will occasionally hand you a malformed record. The adapter absorbs all of
//! that so the tracker never sees a partial contact.

use super::{Batch, Contact, IngestError, SensorHealth};
use crate::geo::{Geodetic, FT_TO_M};
use anyhow::Result;
use serde::Deserialize;
use serde_json::Value;
use std::time::{Duration, Instant};
use tokio::sync::mpsc::Sender;

/// Initial body buffer. Deliberately not the configured limit: reserving the ceiling on
/// every poll would allocate the very memory the ceiling exists to protect.
const BODY_PREALLOC: usize = 64 * 1024;

/// Ceiling on the age a single report is allowed to claim. Anything beyond this is a
/// decode artefact rather than a very patient aircraft, and letting it through would
/// hand the tracker an enormous negative time step to reason about.
const MAX_REPORTED_AGE_S: f64 = 60.0;

#[derive(Debug, Deserialize)]
struct Response {
    #[serde(default)]
    ac: Vec<Aircraft>,
}

#[derive(Debug, Deserialize)]
struct Aircraft {
    hex: Option<String>,
    flight: Option<String>,
    r#type: Option<String>,
    t: Option<String>,
    lat: Option<f64>,
    lon: Option<f64>,
    /// Feet, or the string "ground".
    alt_baro: Option<Value>,
    alt_geom: Option<f64>,
    gs: Option<f64>,
    track: Option<f64>,
    baro_rate: Option<f64>,
    geom_rate: Option<f64>,
    squawk: Option<String>,
    emergency: Option<String>,
    seen_pos: Option<f64>,
}

pub struct AdsbSensor {
    client: reqwest::Client,
    endpoint: String,
    period: Duration,
    /// Hard ceiling on one response body, in bytes.
    max_body: usize,
    /// Owned by this task alone. No `Arc`, no lock — it leaves only by being copied
    /// onto an outgoing batch.
    health: SensorHealth,
}

impl AdsbSensor {
    pub fn new(endpoint: String, period: Duration, max_body: usize) -> Result<Self> {
        let client = reqwest::Client::builder()
            .timeout(Duration::from_secs(10))
            .user_agent("aether/0.1 (telemetry research)")
            .build()?;
        Ok(Self {
            client,
            endpoint,
            period,
            max_body,
            health: SensorHealth::default(),
        })
    }

    /// Poll forever. A failed poll is a dropout, not a fault: it is counted, reported,
    /// and retried on the next tick. A surveillance chain does not get to exit because
    /// one response was unreadable.
    pub async fn run(mut self, tx: Sender<Batch>) {
        let mut tick = tokio::time::interval(self.period);
        tick.set_missed_tick_behavior(tokio::time::MissedTickBehavior::Delay);

        loop {
            tick.tick().await;
            self.health.polls += 1;

            // Either outcome produces a batch. A failed poll sends no contacts, which
            // is not the same statement as "the sky is empty" — the health snapshot
            // riding alongside is what tells the two apart.
            let (received, contacts) = match self.poll().await {
                Ok(pair) => pair,
                Err(e) => {
                    self.health.record(&e);
                    (Instant::now(), Vec::new())
                }
            };

            // A full channel means the tracker is behind, so the batch is dropped
            // rather than blocking the sensor — stale plots are worse than none.
            let _ = tx.try_send(Batch {
                received,
                contacts,
                health: self.health,
            });

            // The only reason to stop polling is that nobody is listening any more.
            if tx.is_closed() {
                return;
            }
        }
    }

    async fn poll(&self) -> Result<(Instant, Vec<Contact>), IngestError> {
        let response = self
            .client
            .get(&self.endpoint)
            // Stated at the boundary rather than inferred from a Cargo feature list.
            // The body limit counts decoded bytes, so an automatically decompressed
            // response would be measured after expansion. This build cannot decompress
            // — `reqwest` is built with `default-features = false` and no gzip, brotli
            // or deflate — but that is a fact about the manifest, and manifests drift.
            // Asking for identity makes the assumption explicit and self-enforcing.
            .header("Accept-Encoding", "identity")
            .send()
            .await
            .map_err(|e| IngestError::Http(e.to_string()))?
            .error_for_status()
            .map_err(|e| IngestError::Http(e.to_string()))?;

        check_declared_length(response.content_length(), self.max_body)?;
        let body = read_bounded(response, self.max_body).await?;

        // Stamped once the body is in hand, which is the closest observable moment to
        // the server's snapshot. Network transit is not corrected for and is left
        // folded into the reported ages — tens of milliseconds against a gate of the
        // order of a hundred metres.
        let received = Instant::now();

        let parsed: Response =
            serde_json::from_slice(&body).map_err(|e| IngestError::Decode(e.to_string()))?;

        Ok((received, parsed.ac.iter().filter_map(decode).collect()))
    }
}

/// Cheap pre-check on the declared length.
///
/// Advisory only, and never the control: `Content-Length` is optional, and a server that
/// wants to flood us can simply omit it or understate it. It buys an early rejection in
/// the honest case, which is the common one.
fn check_declared_length(declared: Option<u64>, limit: usize) -> Result<(), IngestError> {
    match declared {
        Some(len) if len > limit as u64 => Err(IngestError::Oversized {
            limit,
            seen: len.min(usize::MAX as u64) as usize,
        }),
        _ => Ok(()),
    }
}

/// Append one chunk, refusing to let the buffer pass `limit`.
///
/// The size is checked *before* the copy, so the buffer never momentarily holds more
/// than the limit. The bound is on what gets allocated, not on what is noticed
/// afterwards — a check applied after `extend_from_slice` would already have paid for
/// the memory it is supposed to be refusing.
fn accumulate(buf: &mut Vec<u8>, chunk: &[u8], limit: usize) -> Result<(), IngestError> {
    let seen = buf.len() + chunk.len();
    if seen > limit {
        return Err(IngestError::Oversized { limit, seen });
    }
    buf.extend_from_slice(chunk);
    Ok(())
}

/// Read a body chunk by chunk, abandoning it the moment it passes `limit`.
///
/// Returning `Err` here drops the `Response`, which closes the connection: the rest of
/// an oversized body is never read, let alone buffered or parsed. `serde_json` cannot
/// be handed more than `limit` bytes because it is only reached on the `Ok` path.
async fn read_bounded(
    mut response: reqwest::Response,
    limit: usize,
) -> Result<Vec<u8>, IngestError> {
    let mut buf = Vec::with_capacity(BODY_PREALLOC.min(limit));
    while let Some(chunk) = response
        .chunk()
        .await
        .map_err(|e| IngestError::Http(e.to_string()))?
    {
        accumulate(&mut buf, &chunk, limit)?;
    }
    Ok(buf)
}

/// Turn one wire record into a contact, or discard it. A record without an identity or
/// a position is not a degraded contact, it is not a contact.
fn decode(a: &Aircraft) -> Option<Contact> {
    let id = a.hex.as_ref()?.trim().to_ascii_lowercase();
    if id.is_empty() {
        return None;
    }
    let (lat, lon) = (a.lat?, a.lon?);
    if !lat.is_finite() || !lon.is_finite() || lat.abs() > 90.0 || lon.abs() > 180.0 {
        return None;
    }

    let (alt_m, on_ground) = altitude(a)?;

    let label = a
        .flight
        .as_deref()
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .unwrap_or(&id)
        .to_string();

    Some(Contact {
        id,
        label,
        kind: a.t.clone().or_else(|| a.r#type.clone()).unwrap_or_default(),
        squawk: a.squawk.clone().unwrap_or_default(),
        // 7500 hijack, 7600 radio failure, 7700 general emergency.
        emergency: matches!(a.emergency.as_deref(), Some(e) if !e.is_empty() && e != "none")
            || matches!(a.squawk.as_deref(), Some("7500" | "7600" | "7700")),
        on_ground,
        geo: Geodetic {
            lat_deg: lat,
            lon_deg: lon,
            alt_m,
        },
        gs_kt: a.gs.filter(|v| v.is_finite() && *v >= 0.0),
        track_deg: a.track.filter(|v| v.is_finite()),
        vrate_fpm: a.baro_rate.or(a.geom_rate).filter(|v| v.is_finite()),
        age_s: sanitise_age(a.seen_pos),
        source: "ADS-B",
    })
}

/// The reported age of a position, reduced to something the tracker can do arithmetic
/// with. Every time step in the estimator is derived from this number, so a NaN here
/// would silently poison a filter rather than fail loudly.
///
/// Written as an explicit match rather than `clamp`: `f64::clamp` propagates NaN, which
/// is the one input that must not survive. Absent, non-finite, or negative ages are
/// treated as fresh — the innovation gate is the thing that catches a lying sensor, and
/// it can only do that if the arithmetic in front of it stays finite.
fn sanitise_age(raw: Option<f64>) -> f64 {
    match raw {
        Some(v) if v.is_finite() => v.clamp(0.0, MAX_REPORTED_AGE_S),
        _ => 0.0,
    }
}

/// Barometric altitude preferred, geometric as fallback, "ground" as the surface.
fn altitude(a: &Aircraft) -> Option<(f64, bool)> {
    match &a.alt_baro {
        Some(Value::Number(n)) => n.as_f64().map(|ft| (ft * FT_TO_M, false)),
        Some(Value::String(s)) if s.eq_ignore_ascii_case("ground") => Some((0.0, true)),
        _ => a.alt_geom.map(|ft| (ft * FT_TO_M, false)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(json: &str) -> Vec<Contact> {
        let r: Response = serde_json::from_str(json).expect("fixture should parse");
        r.ac.iter().filter_map(decode).collect()
    }

    #[test]
    fn decodes_a_normal_report() {
        let c = parse(
            r#"{"ac":[{"hex":"3cc82c","flight":"VJH357  ","t":"C56X","alt_baro":41000,
                       "gs":464.7,"track":117.56,"baro_rate":0,"squawk":"3714",
                       "emergency":"none","lat":39.365753,"lon":18.74184,"seen_pos":0.451}]}"#,
        );
        assert_eq!(c.len(), 1);
        assert_eq!(c[0].id, "3cc82c");
        assert_eq!(c[0].label, "VJH357");
        assert!(!c[0].emergency);
        assert!(
            (c[0].geo.alt_m - 12_496.8).abs() < 1.0,
            "alt was {}",
            c[0].geo.alt_m
        );
    }

    #[test]
    fn drops_records_without_a_position() {
        let c = parse(r#"{"ac":[{"hex":"abc123","flight":"TEST","alt_baro":10000}]}"#);
        assert!(c.is_empty(), "a contact with no position is not a contact");
    }

    #[test]
    fn handles_ground_altitude_string() {
        let c = parse(r#"{"ac":[{"hex":"deadbe","alt_baro":"ground","lat":37.9,"lon":23.7}]}"#);
        assert_eq!(c.len(), 1);
        assert!(c[0].on_ground);
        assert_eq!(c[0].geo.alt_m, 0.0);
    }

    #[test]
    fn flags_emergency_squawks() {
        let c = parse(
            r#"{"ac":[{"hex":"aaa111","squawk":"7700","lat":37.9,"lon":23.7,"alt_baro":30000}]}"#,
        );
        assert!(c[0].emergency);
    }

    #[test]
    fn rejects_out_of_range_coordinates() {
        let c = parse(r#"{"ac":[{"hex":"bad001","lat":991.0,"lon":23.7,"alt_baro":30000}]}"#);
        assert!(c.is_empty());
    }

    #[test]
    fn reported_age_is_sanitised_before_it_can_reach_the_clock() {
        assert_eq!(sanitise_age(Some(2.5)), 2.5);
        assert_eq!(sanitise_age(None), 0.0);
        assert_eq!(sanitise_age(Some(-3.0)), 0.0);
        assert_eq!(
            sanitise_age(Some(f64::NAN)),
            0.0,
            "NaN must not reach a time step"
        );
        assert_eq!(sanitise_age(Some(f64::INFINITY)), 0.0);
        assert_eq!(sanitise_age(Some(1.0e9)), MAX_REPORTED_AGE_S);
    }

    #[test]
    fn an_empty_feed_is_not_an_error() {
        assert!(parse(r#"{"ac":[]}"#).is_empty());
        assert!(parse(r#"{}"#).is_empty());
    }

    // --- response-body bound ---------------------------------------------------------

    const LIMIT: usize = 4096;

    #[test]
    fn a_body_under_the_limit_is_accepted_and_decoded() {
        let body = br#"{"ac":[{"hex":"3cc82c","flight":"VJH357","alt_baro":41000,
                              "lat":39.36,"lon":18.74,"seen_pos":0.4}]}"#;
        let mut buf = Vec::new();
        assert!(accumulate(&mut buf, body, LIMIT).is_ok());
        assert_eq!(buf.len(), body.len());

        let parsed: Response = serde_json::from_slice(&buf).expect("should parse");
        let contacts: Vec<Contact> = parsed.ac.iter().filter_map(decode).collect();
        assert_eq!(contacts.len(), 1);
        assert_eq!(contacts[0].id, "3cc82c");
    }

    #[test]
    fn a_declared_length_over_the_limit_is_rejected_before_the_body_is_touched() {
        assert_eq!(
            check_declared_length(Some(LIMIT as u64 + 1), LIMIT),
            Err(IngestError::Oversized {
                limit: LIMIT,
                seen: LIMIT + 1
            })
        );
        assert!(check_declared_length(Some(LIMIT as u64), LIMIT).is_ok());
        // Absent is the case that matters: it must fall through to the real control
        // rather than being treated as permission.
        assert!(check_declared_length(None, LIMIT).is_ok());
    }

    #[test]
    fn a_body_of_unknown_length_is_bounded_while_it_is_read() {
        // What a chunked or connection-delimited response looks like: no declared
        // length, arriving a kilobyte at a time, never stopping.
        let chunk = vec![b'x'; 1024];
        let mut buf = Vec::new();
        let mut chunks_taken = 0;
        let err = loop {
            match accumulate(&mut buf, &chunk, LIMIT) {
                Ok(()) => chunks_taken += 1,
                Err(e) => break e,
            }
            assert!(chunks_taken < 100, "the accumulator never refused");
        };
        assert_eq!(
            err,
            IngestError::Oversized {
                limit: LIMIT,
                seen: LIMIT + 1024
            }
        );
        assert_eq!(
            chunks_taken, 4,
            "should stop on the chunk that would overflow"
        );
    }

    #[test]
    fn the_accumulator_never_exceeds_the_bound() {
        // The property that matters: the buffer is bounded at every instant, not merely
        // checked after the fact. A limit enforced after the copy has already paid for
        // the memory it claims to refuse.
        for chunk_size in [1usize, 7, 512, 4095, 4096, 4097, 65_536] {
            let chunk = vec![b'x'; chunk_size];
            let mut buf = Vec::new();
            for _ in 0..64 {
                let _ = accumulate(&mut buf, &chunk, LIMIT);
                assert!(
                    buf.len() <= LIMIT,
                    "buffer reached {} with chunk size {chunk_size}",
                    buf.len()
                );
            }
        }
    }

    #[test]
    fn a_single_oversized_chunk_is_refused_outright() {
        let mut buf = Vec::new();
        let huge = vec![b'x'; LIMIT * 4];
        assert!(accumulate(&mut buf, &huge, LIMIT).is_err());
        assert!(buf.is_empty(), "nothing should have been copied");
    }

    /// Serves an oversized body, with no `Content-Length`, to every connection it gets.
    /// Returns the address to poll. Test-only: the shipped binary has no listener.
    async fn oversized_feed() -> std::net::SocketAddr {
        use tokio::io::{AsyncReadExt, AsyncWriteExt};
        use tokio::net::TcpListener;

        let listener = TcpListener::bind("127.0.0.1:0").await.expect("bind");
        let addr = listener.local_addr().expect("addr");
        tokio::spawn(async move {
            while let Ok((mut sock, _)) = listener.accept().await {
                tokio::spawn(async move {
                    // Consume the request first. Closing a socket with the peer's bytes
                    // still unread sends a reset, and the client sees a transport
                    // failure instead of the oversized body we are trying to serve it.
                    let mut scratch = [0u8; 2048];
                    let _ = sock.read(&mut scratch).await;

                    // No Content-Length: the body is delimited by the close, which is
                    // the case where the declared-length pre-check cannot help and the
                    // chunk accumulator is the only thing standing between the process
                    // and an unbounded allocation.
                    let head = b"HTTP/1.1 200 OK\r\nContent-Type: application/json\r\n\
                                 Connection: close\r\n\r\n";
                    let _ = sock.write_all(head).await;
                    let _ = sock.write_all(&vec![b'x'; 256 * 1024]).await;
                    let _ = sock.flush().await;
                    let _ = sock.shutdown().await;
                });
            }
        });
        addr
    }

    /// Next batch, or fail the test rather than hang the suite.
    async fn next_batch(rx: &mut tokio::sync::mpsc::Receiver<Batch>) -> Batch {
        tokio::time::timeout(Duration::from_secs(10), rx.recv())
            .await
            .expect("sensor went silent")
            .expect("channel closed")
    }

    #[tokio::test]
    async fn an_oversized_body_is_classified_as_oversized_not_as_a_transport_error() {
        // Attribution matters: "the feed is flooding us" and "the feed is unreachable"
        // call for different responses, so they must not collapse into one counter.
        let addr = oversized_feed().await;
        let sensor = AdsbSensor::new(format!("http://{addr}/"), Duration::from_secs(1), LIMIT)
            .expect("sensor");
        match sensor.poll().await {
            Err(IngestError::Oversized { limit, .. }) => assert_eq!(limit, LIMIT),
            other => panic!("expected Oversized, got {other:?}"),
        }
    }

    #[tokio::test]
    async fn an_oversized_response_does_not_terminate_the_polling_loop() {
        let addr = oversized_feed().await;
        let (tx, mut rx) = tokio::sync::mpsc::channel(8);
        let sensor = AdsbSensor::new(format!("http://{addr}/"), Duration::from_millis(50), LIMIT)
            .expect("sensor");
        let task = tokio::spawn(sensor.run(tx));

        // Two consecutive batches: the task is not merely alive, it is still polling.
        let first = next_batch(&mut rx).await;
        let second = next_batch(&mut rx).await;

        assert!(
            first.contacts.is_empty(),
            "a refused response must not fabricate contacts"
        );
        assert!(
            second.health.polls >= 2,
            "loop stopped: {:?}",
            second.health
        );
        assert!(
            second.health.oversized_responses >= 2,
            "faults not attributed to the size limit: {:?}",
            second.health
        );
        assert_eq!(second.health.decode_errors, 0, "parser must never have run");

        task.abort();
    }

    #[test]
    fn health_counts_each_failure_kind_separately() {
        let mut h = SensorHealth::default();
        h.record(&IngestError::Http("refused".into()));
        h.record(&IngestError::Oversized {
            limit: 10,
            seen: 11,
        });
        h.record(&IngestError::Decode("bad json".into()));
        h.record(&IngestError::Oversized {
            limit: 10,
            seen: 99,
        });
        assert_eq!(h.http_errors, 1);
        assert_eq!(h.oversized_responses, 2);
        assert_eq!(h.decode_errors, 1);
        assert_eq!(h.faults(), 4);
    }
}
