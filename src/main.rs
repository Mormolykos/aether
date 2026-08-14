//! Aether — real-time aerospace telemetry engine.
//!
//! Two tasks, one channel, no shared mutable state:
//!
//!   sensor task ──(mpsc, bounded)──▶ tracker task ──▶ operator picture
//!
//! The sensor task owns its socket and the tracker owns every track. Nothing is behind
//! a mutex, so there is no lock for the real-time loop to contend on and no chance of a
//! slow network read stalling the picture. The channel is bounded on purpose: when the
//! tracker falls behind, the sensor sheds the batch instead of growing an unbounded
//! queue of increasingly stale plots.

use aether::config::Config;
use aether::conjunction;
use aether::geo::{Frame, Geodetic};
use aether::ingest::adsb::AdsbSensor;
use aether::ingest::SensorHealth;
use aether::render::{self, Cycles};
use aether::track::TrackStore;
use anyhow::Result;
use std::time::{Duration, Instant};
use tokio::sync::mpsc;

#[tokio::main]
async fn main() -> Result<()> {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| ".env".to_string());
    let cfg = Config::load(&path)?;

    let frame = Frame::new(Geodetic {
        lat_deg: cfg.site_lat,
        lon_deg: cfg.site_lon,
        alt_m: cfg.site_alt_m,
    });

    let endpoint = cfg.adsb_endpoint();
    eprintln!("aether: site {:.4} {:.4}", cfg.site_lat, cfg.site_lon);
    eprintln!("aether: feed {endpoint}");
    eprintln!("aether: acquiring...");

    let (tx, mut rx) = mpsc::channel(8);
    let sensor = AdsbSensor::new(
        endpoint,
        Duration::from_millis(cfg.poll_ms),
        cfg.max_body_bytes,
    )?;
    let sensor_task = tokio::spawn(sensor.run(tx));

    let mut store = TrackStore::new(cfg.process_noise, cfg.meas_var, cfg.gate_sigma);
    let mut cycles = Cycles::default();
    let mut health = SensorHealth::default();
    let timeout = Duration::from_secs(cfg.track_timeout_s);

    let mut tick = tokio::time::interval(Duration::from_millis(cfg.cycle_ms));
    tick.set_missed_tick_behavior(tokio::time::MissedTickBehavior::Delay);

    loop {
        tokio::select! {
            _ = tokio::signal::ctrl_c() => break,
            _ = tick.tick() => {
                let now = Instant::now();
                let t0 = now;

                // Drain everything the sensor has produced since the last cycle. The
                // loop never awaits here: a fixed-rate stage that can block on its
                // input is not a fixed-rate stage.
                let (mut plots, mut polls) = (0usize, 0usize);
                while let Ok(batch) = rx.try_recv() {
                    plots += batch.contacts.len();
                    polls += 1;
                    // Latest snapshot wins: the counters are cumulative, so the newest
                    // batch already contains everything the earlier ones reported.
                    health = batch.health;
                    // The batch carries its own arrival instant. The cycle clock is for
                    // ageing and for the picture; it is not a measurement timestamp.
                    store.ingest(&batch.contacts, &frame, batch.received);
                }

                store.prune(now, timeout);
                let alerts = conjunction::screen(
                    &store, now, cfg.horizon_s, cfg.min_horiz_m, cfg.min_vert_m,
                );

                cycles.record(t0.elapsed(), plots, polls);
                render::draw(
                    &store, &alerts, &cycles, &health, &frame, cfg.display_rows, now,
                );
            }
        }
    }

    sensor_task.abort();
    println!(
        "\naether: {} cycles, {} plots, {} tracks initiated, {} dropped, {} gated, \
         {} superseded, worst cycle {:.2} ms",
        cycles.count,
        cycles.contacts_in,
        store.total_initiated,
        store.total_dropped,
        store.total_gated,
        store.total_superseded,
        cycles.worst.as_secs_f64() * 1000.0,
    );
    println!(
        "aether: sensor {} polls, {} http errors, {} oversized, {} decode errors",
        health.polls, health.http_errors, health.oversized_responses, health.decode_errors,
    );
    Ok(())
}
