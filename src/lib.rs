//! Aether: real-time aerospace telemetry ingestion, tracking and conjunction screening.
//!
//! The pipeline is deliberately the same shape as a ground-based surveillance chain:
//!
//!   sensor adapter -> contact reports -> association -> state estimation
//!                  -> conjunction screening -> operator picture
//!
//! Everything downstream of `ingest` is sensor-agnostic and allocation-stable, so the
//! same core runs against a network feed on a workstation or a serial radar link on an
//! embedded target.

pub mod config;
pub mod conjunction;
pub mod ffi;
pub mod geo;
pub mod ingest;
pub mod kalman;
pub mod render;
pub mod track;
